import argparse
import gzip
import json
import os

from Bio import SwissProt


# Data file is uniprot_sprot_human.dat.gz and uniprot_trembl_human.dat.gz at https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/taxonomic_divisions/.
# We can use SeqIO from Bio to read the file.
# Each record in file will have those attributes: https://biopython.org/docs/1.75/api/Bio.SeqRecord.html
# id, name will be loaded for protein. Ensembl IDs(example: Ensembl:ENST00000372839.7) in dbxrefs will be used to create protein and transcript relationship.





def get_dbxrefs(cross_references):
    dbxrefs = []
    for cross_reference in cross_references:
        database_name = cross_reference[0]
        if database_name == 'EMBL':
            for id in cross_reference[1:3]:
                if id != '-':
                    dbxrefs.append({
                        'name': database_name,
                        'id': id
                    })
        elif database_name in ['RefSeq', 'Ensembl', 'MANE-Select']:
            for item in cross_reference[1:]:
                if item != '-':
                    id = item.split('. ')[0]
                    dbxrefs.append({
                        'name': database_name,
                        'id': id
                    })
        else:
            dbxrefs.append({
                'name': cross_reference[0],
                'id': cross_reference[1]
            })
    dbxrefs.sort(key=lambda x: x['name'])
    return dbxrefs

def get_full_name(description):
    rec_name = None
    description_list = description.split(';')
    for item in description_list:
        if item.startswith('RecName: Full=') or item.startswith('SubName: Full='):
            rec_name = item[14:]
            if ' {' in rec_name:
                rec_name = rec_name[0: rec_name.index(' {')]
            break
    return rec_name

def process_file(organism, source, taxonomy_id, input_file_path, output_file_path):
    with gzip.open(input_file_path, 'rt') as input_file, gzip.open(output_file_path, 'wt') as output_file:
        records = SwissProt.parse(input_file)
        for record in records:
            if record.taxonomy_id == [taxonomy_id]:
                dbxrefs = get_dbxrefs(record.cross_references)
                full_name = get_full_name(record.description)
                to_json = {
                    '_key': record.accessions[0],
                    'name': record.entry_name,
                    'organism': organism,
                    'dbxrefs': dbxrefs,
                    'source': source,
                    'source_url': 'https://www.uniprot.org/help/downloads'
                }
                if full_name:
                    to_json['full_name'] = full_name
                output_file.write(json.dumps(to_json) + '\n')

def main():
    parser = argparse.ArgumentParser(description='Uniprot Parser')
    parser.add_argument('--input-file', type=str, required=True, help='Gzipped .dat file')
    parser.add_argument('--output-file', type=str, required=True, help='Output file path. Output format is jsonl.gz')
    parser.add_argument('--organism', type=str, required=True, choices=['Homo Sapiens', 'Mus Musculus'])
    parser.add_argument('--taxonomy-id', type=str, choices=['9606', '10090'], help='Two taxonomy IDs are allowed: 9606 for Homo sapiens, and 10090 for Mus musculus')
    parser.add_argument('--source', type=str, choices=['UniProtKB/Swiss-Prot', 'UniProtKB/TrEMBL'])
    args = parser.parse_args()
    process_file(args.organism, args.source, args.taxonomy_id, args.input_file, args.output_file)

if __name__ == '__main__':
    main()
