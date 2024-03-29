import os
import gzip
import json
import pickle
import argparse
import multiprocessing

def process_file(file_info):
    file_path, output_dir, dataset, label, source, source_url, entrez_ensembl_dict = file_info
    entrez_id = os.path.basename(file_path)
    ensembl_id = entrez_ensembl_dict.get(entrez_id)

    if ensembl_id:
        output_file = os.path.join(output_dir, f'{entrez_id}.jsonl.gz')
        with open(file_path, 'r') as input_file, gzip.open(output_file, 'wt') as output:
            for line in input_file:
                co_entrez_id, score = line.strip().split()
                co_ensembl_id = entrez_ensembl_dict.get(co_entrez_id)
                if co_ensembl_id:
                    _id = f"{entrez_id}_{co_entrez_id}_{label}"
                    _source = f"genes/{ensembl_id}"
                    _target = f"genes/{co_ensembl_id}"
                    output.write(json.dumps({
                        '_id': _id,
                        '_source': _source,
                        '_target': _target,
                        'label': label,
                        'logit_score': score,
                        'source': source,
                        'source_url': source_url
                    }) + '\n')

def main():
    parser = argparse.ArgumentParser(description='Coxpresdb Parser')
    parser.add-argument('--input-dir', type=str, required=True, help='Input directory')
    parser.add-argument('--output-dir', type=str, required=True, help='Output directory')
    parser.add-argument('--dataset', type=str, required=True, help='Dataset')
    parser.add-argument('--label', type=str, required=True, help='Label')
    parser.add-argument('--source', type=str, required=True, help='Source')
    parser.add-argument('--source-url', type=str, required=True, help='Source URL')
    parser.add-argument('--pickle-file', type=str, required=True, help='Path to entrez-to-ensembl.pkl')
    args = parser.parse_args()

    with open(args.pickle_file, 'rb') as f:
        entrez_ensembl_dict = pickle.load(f)

    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)

    file_paths = [os.path.join(args.input_dir, file) for file in os.listdir(args.input_dir)]
    file_infos = [(file_path, args.output_dir, args.dataset, args.label, args.source, args.source_url, entrez_ensembl_dict)
                  for file_path in file_paths]

    with multiprocessing.Pool(multiprocessing.cpu_count()) as pool:
        pool.map(process_file, file_infos)

if __name__ == '__main__':
    main()
