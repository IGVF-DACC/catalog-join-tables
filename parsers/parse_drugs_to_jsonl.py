import argparse
import csv
import gzip
import json

def parse_drug_row(drug_row):
    xrefs = drug_row[6].split(',')
    drug_ontology_terms = []
    for xref in xrefs:
        if xref.startswith('ChEBI:'):
            drug_ontology_terms.append(
                xref.replace('ChEBI:CHEBI:', 'CHEBI_')
            )
    _key = drug_row[0]
    json_line_token = {
        '_key': _key,
        'name': drug_row[1],
        'drug_ontology_terms': drug_ontology_terms,
        'source': 'pharmGKB',
        'source_url': 'https://www.pharmgkb.org/' + 'chemical/' + _key,
    }
    return json_line_token

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input-file', type=str, required=True, help='Gzipped drugs tsv')
    parser.add_argument('--output-file', type=str, required=True, help='Output filename jsonl.gz format.')
    args = parser.parse_args()
    with gzip.open(args.input_file, 'rt') as infile, gzip.open(args.output_file, 'wt') as outfile:
        drug_csv = csv.reader(infile, delimiter='\t')
        next(drug_csv)
        for drug_row in drug_csv:
            outfile.write(json.dumps(parse_drug_row(drug_row)) + '\n')

if __name__ == '__main__':
    main()
