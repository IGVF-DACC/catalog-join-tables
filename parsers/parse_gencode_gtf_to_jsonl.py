import gzip
import json
import argparse

def parse_info_metadata(info):
    parsed_info = {}
    for key, value in zip(info, info[1:]):
        if key in ['gene_id', 'gene_type', 'gene_name', 'transcript_id', 'transcript_type', 'transcript_name', 'hgnc_id', 'mgi_id']:
            parsed_info[key] = value.replace('"', '').replace(';', '')
    return parsed_info

def get_collection_alias(gene_alias_file_path):
    alias_dict = {}
    with gzip.open(gene_alias_file_path, 'rt') as input:
        next(input)
        for line in input:
            hgnc = ''
            ensembl = ''
            mgi = ''
            (tax_id, gene_id, symbol, locus_tag, synonyms, dbxrefs, chromosome, map_location, description, type_of_gene, symbol_from_nomenclature_authority,
             full_name_from_nomenclature_authority, Nomenclature_status, Other_designations, Modification_date, Feature_type) = line.split('\t')
            split_dbxrefs = dbxrefs.split('|')
            for ref in split_dbxrefs:
                if ref.startswith('HGNC:'):
                    hgnc = ref[5:]
                elif ref.startswith('Ensembl:'):
                    ensembl = ref[8:]
                elif ref.startswith('MGI:'):
                    mgi = ref[4:]
            if ensembl or hgnc or mgi:
                complete_synonyms = []
                complete_synonyms.append(symbol)
                for i in synonyms.split('|'):
                    complete_synonyms.append(i)
                for i in Other_designations.split('|'):
                    complete_synonyms.append(i)
                complete_synonyms.append(symbol_from_nomenclature_authority)
                complete_synonyms.append(full_name_from_nomenclature_authority)
                complete_synonyms = list(set(complete_synonyms))
                if '-' in complete_synonyms:
                    complete_synonyms.remove('-')
                alias = {
                    'alias': complete_synonyms,
                    'entrez': gene_id
                }
                if hgnc:
                    alias.update({'hgnc': hgnc})
                    alias_dict[hgnc] = alias
                if mgi:
                    alias.update({'mgi': mgi})
                    alias_dict[mgi] = alias
                if ensembl:
                    alias_dict[ensembl] = alias
    return alias_dict

def get_hgnc_id(id, info, alias_dict):
    hgnc_id = info.get('hgnc_id')
    if not hgnc_id:
        alias = alias_dict.get(id)
        if alias:
            hgnc_id = alias.get('hgnc')
    return hgnc_id

def get_mgi_id(id, info, alias_dict):
    mgi_id = info.get('mgi_id')
    if not mgi_id:
        alias = alias_dict.get(id)
        if alias:
            mgi_id = alias.get('mgi')
    return mgi_id

def get_alias_by_id(id, hgnc_id, mgi_id, alias_dict):
    for key in [id, hgnc_id, mgi_id]:
        if key in alias_dict:
            return alias_dict[key]
    return None

def parse_gtf_line_to_dictionary(line, alias_dict, version, source_url):
    if line.startswith('#'):
        return None
    split_line = line.strip().split()
    if split_line[2] == 'gene':
        info = parse_info_metadata(split_line[8:])
        gene_id = info['gene_id']
        id = gene_id.split('.')[0]
        hgnc_id = get_hgnc_id(id, info, alias_dict)
        mgi_id = get_mgi_id(id, info, alias_dict)
        alias = get_alias_by_id(id, hgnc_id, mgi_id, alias_dict)
        if gene_id.endswith('_PAR_Y'):
            id = id + '_PAR_Y'
        to_json = {
            '_key': id,
            'gene_id': gene_id,
            'gene_type': info['gene_type'],
            'chr': split_line[0],
            'start:long': int(split_line[3]) - 1,
            'end:long': int(split_line[4]),
            'name': info['gene_name'],
            'source': 'GENCODE',
            'version': version,
            'source_url': source_url
        }
        if hgnc_id:
            to_json.update({'hgnc': hgnc_id})
        if mgi_id:
            to_json.update({'mgi': mgi_id})
        if alias:
            to_json.update({
                'alias': alias['alias'],
                'entrez': alias['entrez']
            })
        return to_json
    return None

def main():
    parser = argparse.ArgumentParser(
        prog='GENCODE GTF Parser',
        description='Parse GENCODE GTF file and generate jsonl output'
    )
    parser.add_argument('-i', '--input_path', required=True,
                        help='input file path')
    parser.add_argument('-o', '--output_path', required=True,
                        help='output file path')
    parser.add_argument('-a', '--gene_alias_file_path', required=True,
                        help='gene alias file path')
    parser.add_argument('-v', '--version', required=True,
                        help='GENCODE version')
    parser.add_argument('-u', '--source_url', required=True,
                        help='GENCODE source URL')
    args = parser.parse_args()

    alias_dict = get_collection_alias(args.gene_alias_file_path)
    with open(args.input_path, 'r') as input_file, gzip.open(args.output_path, 'wb') as output_file:
        for line in input_file:
            parsed_line = parse_gtf_line_to_dictionary(line, alias_dict, args.version, args.source_url)
            if parsed_line:
                output_file.write((json.dumps(parsed_line) + '\n').encode('utf-8'))

if __name__ == '__main__':
    main()
