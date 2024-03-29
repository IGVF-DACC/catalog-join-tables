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

NUMERIC_FIELDS = ['start_position', 'end_position']

FIELDS = [
        'varinfo', 'vid', 'variant_vcf', 'variant_annovar', 'start_position',
        'end_position', 'ref_annovar', 'alt_annovar', 'ref_vcf', 'alt_vcf', 'aloft_value', 'aloft_description',
        'apc_conservation', 'apc_conservation_v2', 'apc_epigenetics_active', 'apc_epigenetics',
        'apc_epigenetics_repressed', 'apc_epigenetics_transcription', 'apc_local_nucleotide_diversity',
        'apc_local_nucleotide_diversity_v2', 'apc_local_nucleotide_diversity_v3', 'apc_mappability', 'apc_micro_rna',
        'apc_mutation_density', 'apc_protein_function', 'apc_protein_function_v2', 'apc_protein_function_v3',
        'apc_proximity_to_coding', 'apc_proximity_to_coding_v2', 'apc_proximity_to_tsstes', 'apc_transcription_factor',
        'bravo_an', 'bravo_af', 'filter_status', 'clnsig', 'clnsigincl', 'clndn', 'clndnincl', 'clnrevstat', 'origin',
        'clndisdb', 'clndisdbincl', 'geneinfo', 'polyphen2_hdiv_score', 'polyphen2_hvar_score', 'mutation_taster_score',
        'mutation_assessor_score', 'metasvm_pred', 'fathmm_xf', 'funseq_value', 'funseq_description',
        'genecode_comprehensive_categoty', 'af_total', 'af_asj_female', 'af_eas_female', 'af_afr_male', 'af_female',
        'af_fin_male', 'af_oth_female', 'af_ami', 'af_oth', 'af_male', 'af_ami_female', 'af_afr', 'af_eas_male', 'af_sas',
        'af_nfe_female', 'af_asj_male', 'af_raw', 'af_oth_male', 'af_nfe_male', 'af_asj', 'af_amr_male', 'af_amr_female',
        'af_amr_sas_female', 'af_fin', 'af_afr_female', 'af_sas_male', 'af_amr', 'af_nfe', 'af_eas', 'af_ami_male',
        'af_fin_female', 'sift_cat', 'sift_val', 'polyphen_cat', 'polyphen_val', 'cadd_rawscore', 'cadd_phred',
        'refseq_category', 'tg_afr', 'tg_all', 'tg_amr', 'tg_eas', 'tg_eur', 'tg_sas'
    ]

def convert_freq_value(value):
        if value == '.':
            value = 0

        try:
            value = float(value)
        except:
            pass

        return value

def parse_metadata(info):
        info_obj = {}
        for pair in info.strip().split(';'):
            try:
                key, value = pair.split('=')
            except:
                if len(pair.split('=')) == 1:
                    key = pair.split('=')[0]
                    value = None

            # example of FREQ value: 'Korea1K:0.9545,0.04545|TOPMED:0.8587|dbGaP_PopFreq:0.9243,0.07566'
            if key == 'FREQ':
                for freq in value.split('|'):
                    freq_name, freq_value = freq.split(':')
                    freq_name = freq_name.lower()
                    values = freq_value.split(',')

                    info_obj[f'freq_{freq_name}_ref'] = convert_freq_value(values[0])

                    if len(values) > 1:
                        info_obj[f'freq_{freq_name}_alt'] = convert_freq_value(values[1])
                    else:
                        if convert_freq_value(values[0]) == 1.0:
                            info_obj[f'freq_{freq_name}_alt'] = 0.0

            # e.g. FAVORFullDB/variant_annovar
            if key.startswith('FAVOR'):
                key = key.split('/')[1].lower()

                if key.lower() not in FIELDS:
                    continue

                if key.startswith('apc') or key.startswith('af') or key.startswith('bravo'):
                    try:
                        value = float(value)
                    except:
                        pass

                if key in NUMERIC_FIELDS:
                    try:
                        value = int(value)
                    except:
                        pass

                info_obj[f'annotation_{key}'] = value

        return info_obj

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
    metadata = parse_metadata(data_line[7])
    parsed_line.update(metadata)
    return parsed_line

def process_file(file_info):
    input_file, input_dir, output_dir, seq_repo_dir = file_info
    input_file_path = os.path.join(input_dir, input_file)
    output_file = input_file.replace('.vcf.gz', '.jsonl.gz')
    output_file_path = os.path.join(output_dir, output_file)

    dp = create_dataproxy(f'seqrepo+file://{seq_repo_dir}')
    seq_repo = SeqRepo(seq_repo_dir)
    translator = Translator(data_proxy=dp)

    with gzip.open(input_file_path, 'rt') as input_file, gzip.open(output_file_path, 'wb') as output_file:
        reader = csv.reader(input_file, delimiter='\t')
        writer = jsonlines.Writer(output_file)
        row = next(reader)
        while not (row[0].startswith('#Chrom') or row[0].startswith('#CHROM')):
            row = next(reader)
        for row in reader:
            writer.write(parse_vcf_line_to_dictionary(row, translator, seq_repo))

def process_directory(input_dir, output_dir, seq_repo_dir):
    files = [f for f in os.listdir(input_dir) if f.endswith('.vcf.gz')]
    file_infos = [(file, input_dir, output_dir, seq_repo_dir) for file in files]
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

    process_directory(args.input_dir, args.output_dir, args.seq_repo_dir)

if __name__ == '__main__':
    main()
