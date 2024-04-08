import csv
import gzip
import json
import hashlib
import argparse

def process_file(input_filepath, output_filepath, organism):
    if organism.lower() == 'human':
        organism = 'Homo sapiens'
    elif organism.lower() == 'mouse':
        organism = 'Mus musculus'
    else:
        raise ValueError("Invalid organism. Please use 'human' or 'mouse'.")

    with gzip.open(input_filepath, 'rt') as interaction_file, gzip.open(output_filepath, 'wt') as output_file:
        interaction_csv = csv.reader(interaction_file)
        next(interaction_csv)  # Skip header row

        for row in interaction_csv:
            if row[3] == 'genetic interference':
                continue

            pmid_url = 'http://pubmed.ncbi.nlm.nih.gov/'
            pmids = [pmid.replace("'", '') for pmid in row[2].replace('[', '').replace(']', '').split(', ')]

            _key = hashlib.sha256('_'.join([row[0], row[1], row[4].replace(':', '_')] + pmids).encode()).hexdigest()
            props = {
                '_key': _key,
                '_from': row[0],
                '_to': row[1],
                'detection_method': row[3],
                'detection_method_code': row[4],
                'interaction_type': row[5],
                'interaction_type_code': row[6],
                'confidence_value_biogrid:long': float(row[7]) if row[7] else None,
                'confidence_value_intact:long': float(row[-2]) if row[-2] else None,
                'source': row[-1],
                'pmids': [pmid_url + pmid for pmid in pmids],
                'organism': organism
            }

            output_file.write(json.dumps(props) + '\n')

def main():
    parser = argparse.ArgumentParser(description='Protein Interaction Parser')
    parser.add_argument('--input-filepath', type=str, required=True, help='Path to the input gzipped CSV file')
    parser.add_argument('--output-filepath', type=str, required=True, help='Path to the output gzipped JSONL file')
    parser.add_argument('--organism', type=str, required=True, choices=['human', 'mouse'], help="Organism ('human' or 'mouse')")
    args = parser.parse_args()

    process_file(args.input_filepath, args.output_filepath, args.organism)

if __name__ == '__main__':
    main()
