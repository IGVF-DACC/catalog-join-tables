import argparse
import json

CHROMOSOMES = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX'] #no chrY or chrM in topLD
LDS_FILE_SUFFIX = '_no_filter_0.2_1000000_LD.csv.gz'
ANNOTATIONS_FILE_SUFFIX = '_no_filter_0.2_1000000_info_annotation.csv.gz'

def make_annotation_file_name(ancestry, chromosome):
    return ancestry + '_' + chromosome + ANNOTATIONS_FILE_SUFFIX

def make_lds_file_name(ancestry, chromosome):
    return ancestry + '_' + chromosome + LDS_FILE_SUFFIX

def complete_file_name_to_s3_path(filename, s3_prefix):
    return s3_prefix + filename

def main(args):
    input_items = []
    for chromosome in CHROMOSOMES:
        input_item = {
                '--annotations_path': complete_file_name_to_s3_path(make_annotation_file_name(args.ancestry, chromosome), args.s3_input_prefix),
                '--lds_path': complete_file_name_to_s3_path(make_lds_file_name(args.ancestry, chromosome), args.s3_input_prefix),
                '--label': 'linkage disequilibrium',
                '--source': 'TopLD',
                '--source_url': 'http://topld.genetics.unc.edu/',
                '--chr': chromosome,
                '--ancestry': args.ancestry,
                '--output_prefix': args.s3_output_prefix,
                '--write_mode': args.write_mode,
        }
        input_items.append(input_item)
    with open(args.output_filepath, 'wt') as fp:
        json.dump(input_items, fp, indent=4)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--s3-input-prefix', type=str)
    parser.add_argument('--s3-output-prefix', type=str)
    parser.add_argument('--ancestry', type=str)
    parser.add_argument('--write-mode', choices=['append', 'overwrite'])
    parser.add_argument('--output-filepath', type=str)
    args = parser.parse_args()
    main(args)

