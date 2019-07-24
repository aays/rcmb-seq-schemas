'''
parental_vcf_filter.py - quick and dirty script to filter
for 'usable' SNPs when calling COs
'''

import argparse
from cyvcf2 import VCF, Writer
from tqdm import tqdm

def args():
    parser = argparse.ArgumentParser(
        description='filter for usable SNPs', 
        usage='python3.5 parental_vcf_filter.py [options]')

    parser.add_argument('-v', '--vcf', required=True,
                        type=str, help='Input VCF')
    parser.add_argument('-g', '--gq_threshold', required=False,
                        type=int, help='GQ threshold [default 30]')
    parser.add_argument('-o', '--out', required=True,
                        type=str, help='File to write to')

    args = parser.parse_args()
    if not args.gq_threshold:
        args.gq_threshold = 30

    return args.vcf, args.gq_threshold, args.out

def parse_bases(bases):
    a, b = bases
    a = list(set(a.split('/')))
    b = list(set(b.split('/')))
    return a, b

def check_usable(record):
    usable = False
    alts = record.ALT
    if len(alts) == 2 and len(set([len(i) for i in alts])) == 1:
        usable = True
    elif len(alts) == 1:
        a, b = parse_bases(record.gt_bases)
        if a != b and len(a) == 1 and len(b) == 1:
            usable = True
    return usable

def parse_sites(vcf, out, gq_threshold):
    total_count, kept_count = 0, 0
    vcf_in = VCF(vcf)
    vcf_out = Writer(out, vcf_in)
    for record in tqdm(vcf_in):
        total_count += 1
        if not all(record.gt_quals > gq_threshold) or not record.is_snp:
            continue
        usable = check_usable(record)
        if usable:
            vcf_out.write_record(record)
            kept_count += 1
        else:
            continue
    print('Filtering complete.')
    print('{kept} of {total} records retained.'.format(kept=kept_count, total=total_count))
    print('New VCF written to {f}'.format(f=out))


def main():
    vcf, gq_threshold, out = args()
    parse_sites(vcf, out, gq_threshold)
    print('Good job!')


if __name__ == '__main__':
    main()

        

