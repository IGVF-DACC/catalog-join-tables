import argparse
import csv
import gzip
import hashlib
import json
import os
import re

from collections import defaultdict
from pathlib import Path

import hgvs.dataproviders.uta
from hgvs.easy import parser
from hgvs.extras.babelfish import Babelfish


VAR_ANNO_FILE_INDEX = {
        'var_drug': {'multiple_drugs': 16, 'alleles': 9, 'comparison_alleles': 20},
        'var_pheno': {'multiple_drugs': 19, 'alleles': 9, 'comparison_alleles': 23},
        'var_fa': {'multiple_drugs': 19, 'alleles': 9, 'comparison_alleles': 21}
    }

SOURCE = 'pharmGKB'

SOURCE_URL_PREFIX = 'https://www.pharmgkb.org/'

def build_variant_id(chr, pos_first_ref_base, ref_seq, alt_seq, assembly='GRCh38'):
    # pos_first_ref_base: 1-based position
    key = '{}_{}_{}_{}_{}'.format(str(chr).replace(
        'chr', '').lower(), pos_first_ref_base, ref_seq, alt_seq, assembly)
    return hashlib.sha256(key.encode()).hexdigest()

def build_variant_id_from_hgvs(hgvs_id, validate=True, assembly='GRCh38'):
    # translate hgvs naming to vcf format e.g. NC_000003.12:g.183917980C>T -> 3_183917980_C_T
    if validate:  # use tools from hgvs, which corrects ref allele if it's wrong
        # got connection timed out error occasionally, could add a retry function
        hdp = hgvs.dataproviders.uta.connect()
        babelfish38 = Babelfish(hdp, assembly_name=assembly)
        try:
            chr, pos_start, ref, alt, type = babelfish38.hgvs_to_vcf(
                parser.parse(hgvs_id))
        except Exception as e:
            print(e)
            return None

        if type == 'sub' or type == 'delins':
            return build_variant_id(chr, pos_start+1, ref[1:], alt[1:])
        else:
            return build_variant_id(chr, pos_start, ref, alt)
        # if no need to validate/query ref allele (e.g. single position substitutions) -> use regex match is quicker
    else:
        if hgvs_id.startswith('NC_'):
            chr = int(hgvs_id.split('.')[0].split('_')[1])
            if chr < 23:
                chr = str(chr)
            elif chr == 23:
                chr = 'X'
            elif chr == 24:
                chr = 'Y'
            else:
                print('Error: unsupported chromosome name.')
                return None

            pos_start = hgvs_id.split('.')[2].split('>')[0][:-1]
            if pos_start.isnumeric():
                ref = hgvs_id.split('.')[2].split('>')[0][-1]
                alt = hgvs_id.split('.')[2].split('>')[1]
                return build_variant_id(chr, pos_start, ref, alt)
            else:
                print('Error: wrong hgvs format.')
                return None
        else:
            print('Error: wrong hgvs format.')
            return None

def load_drug_id_mapping(drug_id_mapping_path):
    # e.g. key: '17-alpha-dihydroequilenin sulfate', value: 'PA166238901'
    drug_id_mapping = {}
    with open(drug_id_mapping_path, 'r') as drug_id_mapfile:
        next(drug_id_mapfile)
        for line in drug_id_mapfile:
            drug_row = line.strip('\n').split('\t')
            drug_id_mapping[drug_row[1]] = drug_row[0]
    return drug_id_mapping

def load_gene_id_mapping(gene_id_mapping_path):
    # e.g. key: 'ABCB1', value: 'ENSG00000085563'
    # a few genes mapped to multiple Ensembl IDs, e.g. SLCO1B3 -> ENSG00000111700, ENSG00000257046
    gene_id_mapping = {}
    with open(gene_id_mapping_path, 'r') as gene_id_mapfile:
        gene_id_csv = csv.reader(gene_id_mapfile, delimiter='\t')
        next(gene_id_csv)
        for gene_id_row in gene_id_csv:
            if gene_id_row[3]:
                gene_id_mapping[gene_id_row[5]] = gene_id_row[3]
    return gene_id_mapping

def load_variant_id_mapping(variant_id_mapping_path):
    # e.g. key: 'rs1000002', value: 'NC_000003.12:g.183917980C>T'
    variant_id_mapping = {}
    with open(variant_id_mapping_path, 'r') as variant_id_mapfile:
        next(variant_id_mapfile)
        for line in variant_id_mapfile:
            variant_row = line.strip('\n').split('\t')
            # i.e. dbSNP id, which is unique for each row
            variant_name = variant_row[1]
            # no ref/alt allele in this column, need to match with ids in synonyms
            if ':' not in variant_row[4]:
                # print('no position info for variant:' + variant_name)
                continue
            position_str = variant_row[4].split(':')[0] + ':g.' + variant_row[4].split(':')[1]
            synonyms = variant_row[-1].split(', ')
            variant_ids = []
            for synonym in synonyms:
                # might miss some del/ins variants
                if synonym.startswith(position_str):
                    if 'del' in synonym:
                        if synonym.split('del')[1]:
                            # rs72552763: NC_000006.12:g.160139851_160139853del, NC_000006.12:g.160139851_160139853delGAT are equivalent
                            synonym = synonym.split('del')[0] + 'del'

                    # no alt allele info in ref ids like NC_000003.12:g.183917980=
                    if '=' not in synonym and synonym not in variant_ids:
                        variant_ids.append(synonym)

            if len(variant_ids) < 1:
                # print(variant_name + ' has no hgvs id.')
                continue
            # elif len(variant_ids) > 1:
                #print(variant_name + ' has multipe hgvs ids.')

            # multiple ids for variants with multiple alternative alleles
            variant_id_mapping[variant_name] = variant_ids
    return variant_id_mapping

def load_study_parameters_mapping(study_parameters_mapping_path):
    study_parameters_mapping = defaultdict(list)  # key: variant annotation ID
    # each variant annotation entry (from one publication though) can have multiple study parameter sets
    with open(study_parameters_mapping_path, 'r') as study_mapfile:
        next(study_mapfile)
        for line in study_mapfile:
            study_row = line.strip('\n').split('\t')
            variant_anno_id = study_row[1]
            p_value = None
            # don't know why split is not working when last few columns are empty...
            if len(study_row) == 17:
                p_value = study_row[11]

            study_parameters_mapping[variant_anno_id].append(
                {
                    'study_parameter_id': study_row[0],
                    'study_type': study_row[2],
                    'study_cases': study_row[3],
                    'study_controls': study_row[4],
                    'p-value': p_value,
                    'biogeographical_groups': study_row[-1]
                }
            )
    return study_parameters_mapping

def match_variant_alleles(variant_hgvs_ids, variant_drug_row, file_prefix):
        # retrieve the studied alleles for variants with multiple alt alleles
        # alleles, comparison_alleles columns have mixed formats accross entries, e.g. # e.g. A; AA + GG; GG; ...
        # assuming each entry only studied single alt allele
        if '>' in variant_hgvs_ids[0]:  # multiple substitutions
            alleles = []
            all_alleles = variant_drug_row[VAR_ANNO_FILE_INDEX[file_prefix]['alleles']] + \
                variant_drug_row[VAR_ANNO_FILE_INDEX[file_prefix]
                                 ['comparison_alleles']]
            for allele in all_alleles:
                if allele in ['A', 'C', 'G', 'T'] and allele not in alleles:
                    alleles.append(allele)

            for variant_hgvs_id in variant_hgvs_ids:
                if variant_hgvs_id.split('>')[1] in alleles:
                    return variant_hgvs_id
        else:  # skip del/ins with multiple alleles for now
            return None

def parse_file(input_path: Path, output_filehandle, drug_id_mapping, gene_id_mapping, variant_id_mapping, study_paramters_mapping):
    filename = input_path.name
    file_prefix = '_'.join(filename.split('_')[:2]) #var_drug, var_anno etc. 
    with gzip.open(input_path, 'rt') as infile:
        variant_drug_csv = csv.reader(infile, delimiter='\t')
        next(variant_drug_csv)
        for variant_drug_row in variant_drug_csv:
            variant_name = variant_drug_row[1]
            # ignore haplotypes like CYP3A4*1, CYP3A4*17, only consider those with rsIDs
            if not variant_name.startswith('rs'):
                continue
            # variant info
            variant_anno_id = variant_drug_row[0]
            variant_name = variant_drug_row[1]
            variant_hgvs_ids = variant_id_mapping.get(variant_name)
            if variant_hgvs_ids is None:
                print(variant_name +
                      ' has no matched variant id.')
                continue
            else:
                if len(variant_hgvs_ids) > 1:
                    # print('multiple allele cases: ' + variant_name + ','.join(variant_hgvs_ids))
                    variant_hgvs_id = match_variant_alleles(variant_hgvs_ids, variant_drug_row, file_prefix)
                    if variant_hgvs_id is None:
                        print('no matched alleles for: ' + variant_name + '\t' + variant_anno_id)
                        continue
                else:
                    variant_hgvs_id = variant_hgvs_ids[0]

                if '>' in variant_hgvs_id:
                    variant_id = build_variant_id_from_hgvs(variant_hgvs_id, False)  # skip validate for simple snvs
                else:
                    variant_id = build_variant_id_from_hgvs(variant_hgvs_id)

                if variant_id is None:
                    print(variant_name + ' failed converting hgvs id.')
                    continue

            # study info
            study_info = study_paramters_mapping.get(variant_anno_id)
            if study_info is None:
                print(variant_anno_id + ' has no matched study info.')
                continue
            # gene info
            # can be multiple genes split by ', ', or empty str for NA cases
            gene_symbols = variant_drug_row[2].split(', ')
            if not variant_drug_row[2]:
                continue

            # drug info
            # had to map drug IDs from drug names, but it is error-prone

            # each row can be associated to multiple drugs, split by ', '
            # while ', ' can also be in part of a single drug name, e.g. PA10390: sulfonamides, urea derivatives

            # the multiple_drugs_flag on column 17 indicates if the association applies to each drug individually or in combination
            # but it is sometimes incorrect
            # -> equals to 'and' when there's a single drug in column 4, or empty when there are multiple drugs

            drug_ids = []
            if not variant_drug_row[3]:
                continue
            # add and/or logic?
            multiple_drugs_flag = variant_drug_row[VAR_ANNO_FILE_INDEX[file_prefix].get('multiple_drugs')]

            if multiple_drugs_flag and ', ' in variant_drug_row[3]:
                # retrieve drug names for rows with double quotes
                # e.g. '"interferon alfa-2a, recombinant", "interferon alfa-2b, recombinant", "ribavirin"'
                drug_names = [d.replace('"', '') for d in re.findall("(\\\".*?\\\")", variant_drug_row[3])]
                # otherwise, no double quotes just split by ', '
                if not drug_names:
                    drug_names = variant_drug_row[3].split(', ')
                for drug_name in drug_names:
                    drug_id = drug_id_mapping.get(drug_name)
                    if drug_id is None:
                        print(drug_name + ' has no matched drug id.')
                    else:
                        drug_ids.append(drug_id)
            else:
                drug_name = variant_drug_row[3]
                drug_id = drug_id_mapping.get(drug_name)
                if drug_id is not None:
                    drug_ids.append(drug_id)
                elif ', ' in drug_name:  # try split the drug names by comma, for cases with likely mis-labeled multiple_drugs_flag
                    drug_names = drug_name.split(', ')
                    for drug_name_split in drug_names:
                        drug_id = drug_id_mapping.get(drug_name_split)
                        if drug_id is None:
                            print(drug_name + ' has no matched drug id.')
                        else:
                            drug_ids.append(drug_id)

            if len(drug_ids) == 0:
                continue
            else:
                for drug_id in drug_ids:
                    edge_key = variant_anno_id + '_' + drug_id
                    for gene_symbol in gene_symbols:
                        gene_id_str = gene_id_mapping.get(gene_symbol)
                        if gene_id_str is None:
                            print(gene_symbol + ' has no matched gene id.')
                        # take care of a few genes mapped to multiple Ensembl IDs
                        # maybe should clear out those cases
                        else:
                            gene_ids = gene_id_str.split(', ')
                            for gene_id in gene_ids:
                                _from = 'variants_drugs/' + edge_key
                                _to = 'genes/' + gene_id
                                second_edge_key = edge_key + '_' + gene_id
                                props = {
                                    'variant_id': variant_id,
                                    '_key': second_edge_key,
                                    '_from': _from,
                                    '_to': _to,
                                    'gene_symbol': gene_symbol,
                                    'source': SOURCE,
                                    'source_url': SOURCE_URL_PREFIX + 'variantAnnotation/' + variant_anno_id
                                }
                                output_filehandle.write(json.dumps(props) + '\n')



def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-d', '--drug-id-mapping-file', type=str, required=True, help='parsing_data_files/pharmGKB_chemicals.tsv')
    parser.add_argument('-g', '--gene-id-mapping-file', type=str, required=True, help='parsing_data_files/pharmGKB_genes.tsv')
    parser.add_argument('-v', '--variant-id-mapping-file', type=str, required=True, help='parsing_data_files/pharmGKB_variants.tsv')
    parser.add_argument('-s', '--study-parameters-mapping-file', type=str, required=True, help='parsing_data_files/pharmGKB_study_parameters.tsv')
    parser.add_argument('-i', '--pharmagkb-input-dir', type=str, required=True, help='directory containing input files with prefix var_')
    parser.add_argument('-o', '--output-file', type=str, required=True, help='output filename, jsonl.gz format')
    args = parser.parse_args()
    drug_id_mapping = load_drug_id_mapping(args.drug_id_mapping_file)
    gene_id_mapping = load_gene_id_mapping(args.gene_id_mapping_file)
    variant_id_mapping = load_variant_id_mapping(args.variant_id_mapping_file)
    study_parameters_mapping = load_study_parameters_mapping(args.study_parameters_mapping_file)
    input_paths = [Path(args.pharmagkb_input_dir + '/' + path) for path in os.listdir(args.pharmagkb_input_dir) if path.startswith('var_')]
    with gzip.open(args.output_file, 'wt') as outfile:
        for path in input_paths:
            parse_file(path, outfile, drug_id_mapping, gene_id_mapping, variant_id_mapping, study_parameters_mapping)


if __name__ == '__main__':
    main()

