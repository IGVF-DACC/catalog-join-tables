import argparse
import gzip
import json
import re

def parse_line(line):
    if line.startswith('#'):
        return None
    split_line = line.split()
    if split_line[10] != 'transcript_id':
        return None
    gene_id = re.sub(r'"|;', '', split_line[9])
    transcript_id = re.sub(r'"|;', '', split_line[11])
    source = gene_id.split('.')[0]
    target = transcript_id.split('.')[0]
    return {
        'gene_id': gene_id,
        'transcript_id': transcript_id,
        'source': source,
        'target': target,
    }

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input-file', type=str, required=True, help='gzipped gtf')
    parser.add_argument('--output-file', type=str, required=True, help='output location, jsonl.gz format')
    args = parser.parse_args()
    with gzip.open(args.input_file, 'rt') as infile, gzip.open(args.output_file, 'wt') as outfile:
        for line in infile:
            parsed_line = parse_line(line)
            if parsed_line:
                outfile.write(json.dumps(parsed_line) + '\n')

if __name__ == '__main__':
    main()
