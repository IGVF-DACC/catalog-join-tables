import gzip
import json
import argparse

def parse_info_metadata(info):
    parsed_info = {}
    for key, value in zip(info, info[1:]):
        if key in ['gene_id', 'gene_type', 'gene_name', 'transcript_id', 'transcript_type', 'transcript_name']:
            parsed_info[key] = value.replace('"', '').replace(';', '')
    return parsed_info

def process_file(input_path, output_path, organism, version, source_url):
    with gzip.open(input_path, 'rt') as input_file, gzip.open(output_path, 'wt') as output_file:
        for line in input_file:
            if line.startswith('#'):
                continue

            data_line = line.strip().split()
            if data_line[2] != 'transcript':
                continue
            data = data_line[:8]
            info = parse_info_metadata(data_line[8:])
            transcript_key = info['transcript_id'].split('.')[0]
            if info['transcript_id'].endswith('_PAR_Y'):
                transcript_key = transcript_key + '_PAR_Y'
            gene_key = info['gene_id'].split('.')[0]
            if info['gene_id'].endswith('_PAR_Y'):
                gene_key = gene_key + '_PAR_Y'
            try:
                if organism == 'HUMAN':
                    label = 'gencode_transcript'
                else:
                    label = 'mm_gencode_transcript'
                props = {
                    'transcript_id': info['transcript_id'],
                    'name': info['transcript_name'],
                    'transcript_type': info['transcript_type'],
                    'chr': data[0],
                    'start': str(int(data[3]) - 1),
                    'end': data[4],
                    'gene_name': info['gene_name'],
                    'source': 'GENCODE',
                    'version': version,
                    'source_url': source_url
                }
                output_file.write(json.dumps({
                    '_key': transcript_key,
                    'label': label,
                    'props': props
                }) + '\n')
            except:
                print(f'fail to process for data: {line}')

def main():
    parser = argparse.ArgumentParser(description='GENCODE GTF Parser')
    parser.add_argument('--input', type=str, required=True, help='Input gzipped GTF file path')
    parser.add_argument('--output', type=str, required=True, help='Output gzipped JSONL file path')
    parser.add_argument('--organism', type=str, default='HUMAN', help='Organism (default: HUMAN)')
    parser.add_argument('--version', type=str, help='GENCODE version (default: v43 for human, vM33 for mouse)')
    parser.add_argument('--source_url', type=str, help='GENCODE source URL (default: based on organism)')
    args = parser.parse_args()

    if args.organism not in ['HUMAN', 'MOUSE']:
        raise ValueError('Invalid organism. Allowed values: HUMAN, MOUSE')

    if args.version is None:
        args.version = 'v43' if args.organism == 'HUMAN' else 'vM33'

    if args.source_url is None:
        args.source_url = 'https://www.gencodegenes.org/human/' if args.organism == 'HUMAN' else 'https://www.gencodegenes.org/mouse/'

    process_file(args.input, args.output, args.organism, args.version, args.source_url)

if __name__ == '__main__':
    main()
