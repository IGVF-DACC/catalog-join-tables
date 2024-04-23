import argparse
import csv
import gzip
import hashlib
import json
import logging
import os
import multiprocessing

from math import floor
from math import isinf
from math import log10


SOURCE = 'GTEx'
EQTL_SOURCE_URL_PREFIX = 'https://storage.googleapis.com/adult-gtex/bulk-qtl/v8/single-tissue-cis-qtl/GTEx_Analysis_v8_eQTL/'
SQTL_SOURCE_URL_PREFIX = 'https://storage.googleapis.com/adult-gtex/bulk-qtl/v8/single-tissue-cis-qtl/GTEx_Analysis_v8_sQTL/'
MAX_LOG10_PVALUE = 400

def load_ontology_mapping(ontology_id_mapping_path):
        ontology_id_mapping = {}  # e.g. key: 'Brain_Amygdala', value: 'UBERON_0001876'
        # e.g. filename: Esophagus_Gastroesophageal_Junction -> tissue name: Esophagus - Gastroesophageal Junction
        tissue_name_mapping = {}
        # e.g. key: 'Brain_Amygdala', value (UBERON term name): 'amygdala'
        ontology_term_mapping = {}

        with open(ontology_id_mapping_path, 'r') as ontology_id_mapfile:
            ontology_id_csv = csv.reader(ontology_id_mapfile, delimiter='\t')
            next(ontology_id_csv)
            for row in ontology_id_csv:
                if row[1]:
                    ontology_id_mapping[row[1]] = row[2].replace(':', '_')
                    tissue_name_mapping[row[1]] = row[0]
                    ontology_term_mapping[row[1]] = row[3]
        return ontology_id_mapping, tissue_name_mapping, ontology_term_mapping

def build_variant_id(chr, pos_first_ref_base, ref_seq, alt_seq, assembly='GRCh38'):
    # pos_first_ref_base: 1-based position
    key = '{}_{}_{}_{}_{}'.format(str(chr).replace(
        'chr', '').lower(), pos_first_ref_base, ref_seq, alt_seq, assembly)
    return hashlib.sha256(key.encode()).hexdigest()

def to_float(str):
    MAX_EXPONENT = 307
    number = float(str)

    if number == 0:
        return number
    if isinf(number) and number > 0:
        return float('1e307')
    if isinf(number) and number < 0:
        return float('1e-307')

    base10 = log10(abs(number))
    exponent = floor(base10)

    if abs(exponent) > MAX_EXPONENT:
        if exponent < 0:
            number = number * float(f'1e{abs(exponent) - MAX_EXPONENT}')
        else:
            number = number / float(f'1e{abs(exponent) - MAX_EXPONENT}')
    return number

def parse_eqtl_row_to_jsonl(row, biological_context, input_file):
    chrom, pos, ref_seq, alt_seq, assembly_code = row[0].split('_')
    if assembly_code != 'b38':
        logging.info(f'Unsupported assembly {assembly_code}')
        return None
    variant_id = build_variant_id(chrom, pos, ref_seq, alt_seq)
    gene_id = row[1].split('.')[0]
    variants_genes_id = hashlib.sha256((variant_id + '_' + gene_id + '_' + biological_context).encode()).hexdigest()
    _id = variants_genes_id
    _source = variant_id
    _target = gene_id
    source_url = EQTL_SOURCE_URL_PREFIX + input_file
    pvalue = float(row[6])
    if pvalue == 0:
        log_pvalue = MAX_LOG10_PVALUE
    else:
        log_pvalue = -1 * log10(pvalue)
    jsonl_token = {
        '_id': _id,
        '_source': _source,
        '_target': _target,
        'biological_context': biological_context,
        'log10pvalue': log_pvalue,
        'effect_size': to_float(row[7]),
        'pval_beta': to_float(row[-1]),
        'label': 'eQTL',
        'source': SOURCE,
        'source_url': source_url
    }
    return jsonl_token

def parse_sqtl_row_to_jsonl(row, biological_context, input_file):
    variant_id_info = row[0]
    variant_id_ls = variant_id_info.split('_')
    chrom, pos, ref_seq, alt_seq, _ = variant_id_ls
    variant_id = build_variant_id(chrom, pos, ref_seq, alt_seq)
    phenotype_id = row[1]
    phenotype_id_ls = phenotype_id.split(':')
    gene_id = phenotype_id_ls[-1].split('.')[0]
    variants_genes_id = hashlib.sha256((variant_id + '_' + phenotype_id + '_' + biological_context).encode()).hexdigest()
    _id = variants_genes_id
    _source = variant_id
    _target = gene_id
    source_url = SQTL_SOURCE_URL_PREFIX + input_file
    pvalue = float(row[6])

    if pvalue == 0:
        log_pvalue = MAX_LOG10_PVALUE
    else:
        log_pvalue = -1 * log10(pvalue)

    jsonl_token = {
        '_id': _id,
        '_source': _source,
        '_target': _target,
        'biological_context': biological_context,
        'chr': chrom,
        'sqrt_maf': to_float(row[5]),
        'log10pvalue': log_pvalue,
        'pval_nominal_threshold': to_float(row[9]),
        'min_pval_nominal': to_float(row[10]),
        'p_value': pvalue,
        'effect_size': to_float(row[7]),
        'effect_size_se': to_float(row[9]),
        'pval_beta': to_float(row[11]),
        'intron_chr': phenotype_id_ls[0],
        'intron_start': phenotype_id_ls[1],
        'intron_end': phenotype_id_ls[2],
        'label': 'splice_QTL',
        'source': SOURCE,
        'source_url': source_url,
    }
    return jsonl_token




def process_file(file_info):
    input_file, input_dir, output_dir, ontology_mapping_file, qtl_type = file_info
    _, _ , ontology_term_mapping = load_ontology_mapping(ontology_mapping_file)
    input_file_path = os.path.join(input_dir, input_file)
    output_file = input_file.replace('.txt.gz', '.jsonl.gz')
    output_file_path = os.path.join(output_dir, output_file)
    filename_biological_context = input_file.split('.')[0]
    biological_context = ontology_term_mapping.get(filename_biological_context)
    with gzip.open(input_file_path, 'rt') as infile, gzip.open(output_file_path, 'wt') as outfile:
        reader = csv.reader(infile, delimiter='\t')
        next(reader)
        if qtl_type == 'eqtl':
            row_parser = parse_eqtl_row_to_jsonl
        elif qtl_type == 'sqtl':
            row_parser = parse_sqtl_row_to_jsonl
        else:
            raise ValueError('qtl-type must be eqtl or sqtl')
        for row in reader:
            outfile.write(json.dumps(row_parser(row, biological_context, input_file)) + '\n')

def process_directory(input_dir, output_dir, ontology_mapping_file, qtl_type):
    if qtl_type == 'eqtl':
        file_suffix = 'signif_variant_gene_pairs.txt.gz'
    elif qtl_type == 'sqtl':
        file_suffix = 'sqtl_signifpairs.txt.gz'
    else:
        raise ValueError('qtl-type must be eqtl or sqtl')
    files = [f for f in os.listdir(input_dir) if f.endswith(file_suffix)]
    file_infos = [(file, input_dir, output_dir, ontology_mapping_file, qtl_type) for file in files]
    with multiprocessing.Pool(multiprocessing.cpu_count()) as pool:
        pool.map(process_file, file_infos)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--qtl-type', type=str, required=True, help='eqtl or sqtl')
    parser.add_argument('--ontology-mapping-file', type=str, required=True, help='Ontology info, parsing_data_files/GTEx_UBERON_mapping.tsv')
    parser.add_argument('--input-dir', type=str, required=True, help='Suffixes sqtl_signifpairs.txt.gz in case of sqtl, and signif_variant_gene_pairs.txt.gz in case of eqtl will be processed')
    parser.add_argument('--output-dir', type=str, required=True, help='One jsonl.gz output will be written per one input file into this directory')
    args = parser.parse_args()
    process_directory(args.input_dir, args.output_dir, args.ontology_mapping_file, args.qtl_type)

if __name__ == '__main__':
    main()
