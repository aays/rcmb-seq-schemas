'''
read_snp_counts.py - get snp counts for different read-pair schemes
'''

import argparse
from tqdm import tqdm
import ah_utils
import csv
from cyvcf2 import VCF

def args():
    parser = argparse.ArgumentParser(
        description='get SNP counts for different read-pair schemes', 
        usage='python3.5 read_snp_counts.py [options]')

    parser.add_argument('-v', '--vcf', required=True,
                        type=str, help='Filtered input VCF')
    parser.add_argument('-l', '--read_length', required=True,
                        type=int, help='Read length')
    parser.add_argument('-i', '--insert_size', required=True,
                        type=int, help='Insert size b/w reads')
    parser.add_argument('-c', '--chrom', required=True,
                        type=str, help='Chromosome')
    parser.add_argument('-o', '--out', required=True,
                        type=str, help='File to write to')

    args = parser.parse_args()

    return args.vcf, args.read_length, args.insert_size, \
           args.chrom, args.out

def count_window_snps(vcf, start, read_length, insert_size, chrom):
    ''' (str, int, int, int, str) -> int, int
    helper function for parse_windows below
    counts SNPs in both left and right read given start point + read_length + ins size
    '''
    window_left = (start, start + read_length)
    right_start = window_left[1] + insert_size
    window_right = (right_start, right_start + read_length)
    left_count, right_count = 0, 0

    vcfin = VCF(vcf)

    range_str = '{0}:{1}-{2}'
    for record in vcfin.__call__(range_str.format(chrom, *window_left)):
        left_count += 1
    for record in vcfin.__call__(range_str.format(chrom, *window_right)):
        right_count += 1
    
    return left_count, right_count


def parse_windows(vcf, read_length, insert_size, chrom, out):
    ''' (str, int, int, str, str) -> None
    iterates through VCF, providing counts of SNPs at each read for
    given read_length + insert_size pairing
    '''

    chrom_length = ah_utils.chlamy_lengths()[chrom]

    with open(out, 'w') as f:
        fieldnames = [
            'chrom', 'read_length' ,'insert_size', 'left_start', 
            'left_end', 'right_start', 'right_end', 'left_count', 'right_count'
            
        ]
        writer = csv.DictWriter(f, delimiter='\t', fieldnames=fieldnames)
        writer.writeheader()

        windows = list(range(0, chrom_length, int(read_length / 2)))
        windows.append(chrom_length)

        for start in tqdm(windows):
            left_count, right_count = count_window_snps(
                vcf, start, read_length, insert_size, chrom
            )
            out_dict = dict.fromkeys(fieldnames)
            out_dict['chrom'] = chrom
            out_dict['read_length'] = read_length
            out_dict['insert_size'] = insert_size
            out_dict['left_start'] = start
            out_dict['left_end'] = start + read_length
            out_dict['right_start'] = start + read_length + insert_size
            out_dict['right_end'] = start + insert_size + (2 * read_length)
            out_dict['left_count'] = left_count
            out_dict['right_count'] = right_count

            writer.writerow(out_dict)


def main():
    vcf, read_length, insert_size, chrom, out = args()
    print('SNP counts for L {0} with I {1}'.format(read_length, insert_size))
    parse_windows(vcf, read_length, insert_size, chrom, out)
    

if __name__ == '__main__':
    main()

        

