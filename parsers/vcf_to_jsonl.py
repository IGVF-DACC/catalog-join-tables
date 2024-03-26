import csv
import hashlib
import argparse
import gzip
import os
import multiprocessing

import jsonlines

from ga4gh.vrs.extras.translator import Translator, ValidationError
from ga4gh.vrs.dataproxy import create_dataproxy
from biocommons.seqrepo import SeqRepo

CHR_MAP = {
    '1': 'NC_000001.11',
    '2': 'NC_000002.12',
    '3': 'NC_000003.12',
    '4': 'NC_000004.12',
    '5': 'NC_000005.10',
    '6': 'NC_000006.12',
    '7': 'NC_000007.14',
    '8': 'NC_000008.11',
    '9': 'NC_000009.12',
    '10': 'NC_000010.11',
    '11': 'NC_000011.10',
    '12': 'NC_000012.12',
    '13': 'NC_000013.11',
    '14': 'NC_000014.9',
    '15': 'NC_000015.10',
    '16': 'NC_000016.10',
    '17': 'NC_000017.11',
    '18': 'NC_000018.10',
    '19': 'NC_000019.10',
    '20': 'NC_000020.11',
    '21': 'NC_000021.9',
    '22': 'NC_000022.11',
    'X': 'NC_000023.11',
    'Y': 'NC_000024.10'
}

def build_variant_id(chr, pos_first_ref_base, ref_seq, alt_seq, assembly='GRCh38'):
    key = '{}_{}_{}_{}_{}'.format(str(chr).replace(
        'chr', '').lower(), pos_first_ref_base, ref_seq, alt_seq, assembly)
    return hashlib.sha256(key.encode()).hexdigest()

def build_allele(chr, pos, ref, alt, translator, seq_repo):
    gnomad_exp = f'{chr}-{pos}-{ref}-{alt}'
    try:
        allele = translator.translate_from(gnomad_exp, 'gnomad')
    except ValidationError as e:
        print(e)
        chr_ref = CHR_MAP[chr]
        start = int(pos) - 1
        end = start + len(ref)
        ref = seq_repo[chr_ref][start:end]
        gnomad_exp = f'{chr}-{pos}-{ref}-{alt}'
        print('correct gnomad_exp:', gnomad_exp)
        allele = translator.translate_from(gnomad_exp, 'gnomad')
    return allele

def build_spdi(chr, pos, ref, alt, translator, seq_repo):
    if len(ref) == 1 and len(alt) == 1:
        chr_ref = CHR_MAP[chr]
        pos_spdi = int(pos) - 1
        spdi = f'{chr_ref}:{pos_spdi}:{ref}:{alt}'
    else:
        allele = build_allele(chr, pos, ref, alt, translator, seq_repo)
        spdi = translator.translate_to(allele, 'spdi')[0]
        del_seq = translator.data_proxy.get_sequence(str(
            allele.location.sequence_id), allele.location.interval.start.value, allele.location.interval.end.value)
        spdi = convert_spdi(spdi, del_seq)
    return spdi

def convert_spdi(spdi, seq):
    ls = spdi.split(':')
    ls[2] = seq
    spdi = ':'.join(ls)
    return spdi

def build_hgvs_from_spdi(spdi):
    ins_seq = spdi.split(':')[-1]
    del_seq = spdi.split(':')[2]
    spdi_pos = int(spdi.split(':')[1])
    chr_ref = spdi.split(':')[0]
    if len(ins_seq) == 1 and len(del_seq) == 1:
        hgvs = f'{chr_ref}:g.{spdi_pos + 1}{del_seq}>{ins_seq}'
    elif del_seq.startswith(ins_seq):
        pos_hgvs_start = spdi_pos + 1 + len(ins_seq)
        pos_hgvs_end = spdi_pos + len(del_seq)
        if pos_hgvs_start == pos_hgvs_end:
            hgvs = f'{chr_ref}:g.{pos_hgvs_start}del'
        else:
            hgvs = f'{chr_ref}:g.{pos_hgvs_start}_{pos_hgvs_end}del'
    elif ins_seq.startswith(del_seq):
        pos_hgvs_start = spdi_pos + len(del_seq)
        pos_hgvs_end = pos_hgvs_start + 1
        insert_seq_hgvs = ins_seq[len(del_seq):]
        hgvs = f'{chr_ref}:g.{pos_hgvs_start}_{pos_hgvs_end}ins{insert_seq_hgvs}'
    else:
        pos_hgvs_start = spdi_pos + 1
        pos_hgvs_end = spdi_pos + len(del_seq)
        if pos_hgvs_start == pos_hgvs_end:
            hgvs = hgvs = f'{chr_ref}:g.{pos_hgvs_start}delins{ins_seq}'
        else:
            hgvs = f'{chr_ref}:g.{pos_hgvs_start}_{pos_hgvs_end}delins{ins_seq}'
    return hgvs

def parse_vcf_line_to_dictionary(data_line, translator, seq_repo) -> dict:
    spdi = build_spdi(data_line[0], data_line[1], data_line[3], data_line[4], translator, seq_repo)
    hgvs = build_hgvs_from_spdi(spdi)
    parsed_line = {
        'chr': f'chr{data_line[0]}',
        'pos:long': int(data_line[1]) - 1,
        'rsid': data_line[2],
        'ref': data_line[3],
        'alt': data_line[4],
        'spdi': spdi,
        'hgvs': hgvs,
        'source': 'FAVOR',
        'source_url': 'http://favor.genohub.org/',
    }
    return parsed_line

def process_file(file_info):
    input_file, input_dir, output_dir, translator, seq_repo = file_info
    input_file_path = os.path.join(input_dir, input_file)
    output_file = input_file.replace('.vcf.gz', '.jsonl.gz')
    output_file_path = os.path.join(output_dir, output_file)
    with gzip.open(input_file_path, 'rt') as input_file, gzip.open(output_file_path, 'wb') as output_file:
        reader = csv.reader(input_file, delimiter='\t')
        writer = jsonlines.Writer(output_file)
        row = next(reader)
        while not (row[0].startswith('#Chrom') or row[0].startswith('#CHROM')):
            row = next(reader)
        for row in reader:
            writer.write(parse_vcf_line_to_dictionary(row, translator, seq_repo))

def process_directory(input_dir, output_dir, translator, seq_repo):
    files = [f for f in os.listdir(input_dir) if f.endswith('.vcf.gz')]
    file_infos = [(file, input_dir, output_dir, translator, seq_repo) for file in files]
    with multiprocessing.Pool(multiprocessing.cpu_count()) as pool:
        pool.map(process_file, file_infos)

def main():
    parser = argparse.ArgumentParser(
        prog='Variants SPDI generator',
        description='Generate SPDI for variants'
    )
    parser.add_argument('-i', '--input_dir', required=True,
                        help='input directory path')
    parser.add_argument('-o', '--output_dir', required=True,
                        help='output directory path')
    parser.add_argument('--seq-repo-dir', required=True, help='Seqrepo location')
    args = parser.parse_args()

    dp = create_dataproxy(f'seqrepo+file://{args.seq_repo_dir}')
    seq_repo = SeqRepo(args.seq_repo_dir)
    translator = Translator(data_proxy=dp)

    process_directory(args.input_dir, args.output_dir, translator, seq_repo)

if __name__ == '__main__':
    main()

