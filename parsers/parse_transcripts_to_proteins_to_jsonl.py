import argparse
import gzip
import json

from Bio import SeqIO

'''
Data file is uniprot_sprot_human.dat.gz and uniprot_trembl_human.dat.gz at https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/taxonomic_divisions/.
We can use SeqIO from Bio to read the file.
Each record in file will have those attributes: https://biopython.org/docs/1.75/api/Bio.SeqRecord.html
id, name will be loaded for protein. Ensembl IDs(example: Ensembl:ENST00000372839.7) in dbxrefs will be used to create protein and transcript relationship.
'''

ORGANISM_TO_ENSEMBL_PREFIX = {
    'human': 'ENST',
    'mouse': 'ENSMUST',
}


def process_record(record, ensembl_prefix):
    dbxrefs = record.dbxrefs
    results = []
    for item in dbxrefs:
        if item.startswith('Ensembl') and ensembl_prefix in item:
            ensg_id = item.split(':')[-1].split('.')[0]
            _id = ensg_id + '_' + record.id
            _source = ensg_id
            _target = record.id
            token = {
                'id': _id,
                'source': _source,
                'target':  _target,
            }
            results.append(token)
    return results


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--organism', type=str, required=True, default='human', help='human or mouse, default: human')
    parser.add_argument('--input-file', type=str, required=True, help='Input .dat.gz')
    parser.add_argument('--output-file', type=str, required=True, help='Output file, .jsonl.gz')
    args = parser.parse_args()
    if args.organism.lower() not in ['human', 'mouse']:
        raise ValueError('--organism must be human or mouse')

    ensembl_prefix = ORGANISM_TO_ENSEMBL_PREFIX[args.organism.lower()]

    with gzip.open(args.input_file, 'rt') as infile, gzip.open(args.output_file, 'wt') as outfile:
        records = SeqIO.parse(infile, 'swiss')
        for record in records:
            if not record.name.endswith(args.organism.upper()):
                continue
            else:
                jsonl_tokens_from_record = process_record(record, ensembl_prefix)
                for token in jsonl_tokens_from_record:
                    outfile.write(json.dumps(token) + '\n')

if __name__ == '__main__':
    main()
