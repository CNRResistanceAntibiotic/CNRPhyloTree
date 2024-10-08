#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
    This script will manage the vcf file to matrix distance program
"""

import os
import argparse
from collections import OrderedDict
import vcf
from itertools import combinations
import pandas as pd


def read_vcf(vcf_file):
    print(f"\nRead vcf file {vcf_file}")
    vcf_dic = OrderedDict()
    vcf_reader = vcf.Reader(open(vcf_file, 'r'))
    sample_names = vcf_reader.samples
    for record in vcf_reader:
        if record.CHROM not in vcf_dic.keys():
            vcf_dic[record.CHROM] = [record]
        else:
            vcf_dic[record.CHROM].append(record)
    print(f"Number of loci   : {len(vcf_dic)}")
    print(f"Number of samples: {len(sample_names)}")
    return vcf_dic, sample_names


def vcf2dist(vcf_dic, sample_names, type_matrice):
    print("\nConvert VCF data in SNP and MLST distances")
    contig_dist_dict = OrderedDict()
    snp_dist_dict = OrderedDict()
    less_snp_dist_dict = OrderedDict()
    count_record = 0
    for sample1 in sample_names:
        if sample1 not in contig_dist_dict.keys():
            contig_dist_dict[sample1] = OrderedDict()
            snp_dist_dict[sample1] = OrderedDict()
            less_snp_dist_dict[sample1] = OrderedDict()
        for sample2 in sample_names:
            contig_dist_dict[sample1][sample2] = 0
            snp_dist_dict[sample1][sample2] = 0
            less_snp_dist_dict[sample1][sample2] = 0

    for chrom in vcf_dic.keys():
        records = vcf_dic[chrom]
        mlst = OrderedDict()
        for sample1, sample2 in combinations(sample_names, 2):
            if sample1 not in mlst.keys():
                mlst[sample1] = OrderedDict({sample2: 0})
            else:
                mlst[sample1][sample2] = 0
        # snp distance
        for record in records:
            count_record += 1
            for sample1, sample2 in combinations(sample_names, 2):
                # print sample1, record.genotype(sample1)['GT'], sample2, record.genotype(sample2)['GT']
                if record.genotype(sample1)['GT'] != record.genotype(sample2)['GT']:
                    s1_genotype = str(record.ALT[int(record.genotype(sample1)['GT'].split('/')[0]) - 1])
                    s2_genotype = str(record.ALT[int(record.genotype(sample2)['GT'].split('/')[0]) - 1])

                    if s1_genotype != 'N' and s2_genotype != 'N':
                        snp_dist_dict[sample1][sample2] = snp_dist_dict[sample1][sample2] + 1
                        snp_dist_dict[sample2][sample1] = snp_dist_dict[sample1][sample2]
                        mlst[sample1][sample2] = mlst[sample1][sample2] + 1
                    else:
                        less_snp_dist_dict[sample1][sample2] = less_snp_dist_dict[sample1][sample2] + 1
                        less_snp_dist_dict[sample2][sample1] = less_snp_dist_dict[sample1][sample2]
                        print(f"{sample1} {sample2} {chrom} {record.POS}: N genotype excluded in SNP count")
        # contig distance (unscoop to snp distance)
        for sample1, sample2 in combinations(sample_names, 2):
            if mlst[sample1][sample2] > 0:
                contig_dist_dict[sample1][sample2] = contig_dist_dict[sample1][sample2] + 1
                contig_dist_dict[sample2][sample1] = contig_dist_dict[sample1][sample2]

    if type_matrice == "relatif-norm":
        for sample1, sample2 in combinations(sample_names, 2):
            value_to_subtract = less_snp_dist_dict[sample1][sample2]
            denominator = count_record-value_to_subtract
            value_init = snp_dist_dict[sample1][sample2]
            snp_dist_dict[sample1][sample2] = (value_init/denominator) * 100
            snp_dist_dict[sample2][sample1] = (value_init / denominator) * 100

    return snp_dist_dict, contig_dist_dict


def write_dist(dist_dic, dist_type, out_prefix, sample_names):
    df = pd.DataFrame.from_records(data=dist_dic, columns=sample_names)
    outfile = out_prefix + f'_{dist_type}_dist.tsv'
    df.to_csv(outfile, sep='\t')
    print(f"{dist_type} distance matrix written in {outfile}")


def pre_main(args):
    vcf_file = os.path.abspath(args.vcf)
    type_matrice = args.type_matrice
    main(vcf_file, type_matrice)


def main(vcf_file, type_matrice):
    out_prefix = os.path.splitext(vcf_file)[0]
    vcf_dic, sample_names = read_vcf(vcf_file)
    snp_dist_dict, contig_dist_dict = vcf2dist(vcf_dic, sample_names, type_matrice)
    write_dist(snp_dist_dict, 'SNP', out_prefix, sample_names)
    write_dist(contig_dist_dict, 'Contig', out_prefix, sample_names)
    print('\nvcf to matrix distances done!\n')


def version():
    return '1.0'


def run():
    parser = argparse.ArgumentParser(description='vcf2dist - Version ' + version())
    parser.add_argument('-v', '--vcf', dest="vcf", default='output_filtered.vcf', help='vcf file [output_filtered.vcf]')
    parser.add_argument('-V', '--version', action='version', version='vcf2dist-' + version(),
                        help="Prints version number")
    parser.add_argument('-tm', '--type_matrice', dest="type_matrice", default="absolu",
                        help="The mode of computation of the distance matrice. [absolu, relatif, relatif-norm]"
                             " (default: absolu)")
    args = parser.parse_args()
    pre_main(args)


if __name__ == '__main__':
    run()
