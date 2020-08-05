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
    print("\nRead vcf file {}".format(vcf_file))
    vcf_dic = OrderedDict()
    vcf_reader = vcf.Reader(open(vcf_file, 'r'))
    sample_names = vcf_reader.samples
    for record in vcf_reader:
        if record.CHROM not in vcf_dic.keys():
            vcf_dic[record.CHROM] = [record]
        else:
            vcf_dic[record.CHROM].append(record)
    print("Number of loci   : {}".format(len(vcf_dic)))
    print("Number of samples: {}".format(len(sample_names)))
    return vcf_dic, sample_names


def vcf2dist(vcf_dic, sample_names):
    print("\nConvert VCF data in SNP and MLST distances")
    contig_dist_dict = OrderedDict()
    snp_dist_dict = OrderedDict()
    contig_dist_rework_dict = OrderedDict()
    snp_dist_rework_dict = OrderedDict()

    for sample1 in sample_names:
        if sample1 not in contig_dist_dict.keys():
            contig_dist_dict[sample1] = OrderedDict()
            snp_dist_dict[sample1] = OrderedDict()

        if sample1 not in contig_dist_rework_dict.keys():
            contig_dist_rework_dict[sample1] = OrderedDict()
            snp_dist_rework_dict[sample1] = OrderedDict()

        for sample2 in sample_names:
            contig_dist_dict[sample1][sample2] = 0
            snp_dist_dict[sample1][sample2] = 0

            contig_dist_rework_dict[sample1][sample2] = 0
            snp_dist_rework_dict[sample1][sample2] = 0

    for chrom in vcf_dic.keys():
        records = vcf_dic[chrom]
        mlst = OrderedDict()
        mlst_rework = OrderedDict()
        for sample1, sample2 in combinations(sample_names, 2):
            if sample1 not in mlst.keys():
                mlst[sample1] = OrderedDict({sample2: 0})
            else:
                mlst[sample1][sample2] = 0

            if sample1 not in mlst_rework.keys():
                mlst_rework[sample1] = OrderedDict({sample2: 0})
            else:
                mlst_rework[sample1][sample2] = 0

        # snp distance
        for record in records:
            for sample1, sample2 in combinations(sample_names, 2):
                # print sample1, record.genotype(sample1)['GT'], sample2, record.genotype(sample2)['GT']
                if record.genotype(sample1)['GT'] != record.genotype(sample2)['GT']:
                    # print 'sample1 diff sample2'
                    s1_genotype = str(record.ALT[int(record.genotype(sample1)['GT'].split('/')[0]) - 1])
                    s2_genotype = str(record.ALT[int(record.genotype(sample2)['GT'].split('/')[0]) - 1])
                    # print s1_genotype, s2_genotype
                    if s1_genotype != 'N' or s2_genotype != 'N':
                        # print 'diff of N'
                        snp_dist_dict[sample1][sample2] = snp_dist_dict[sample1][sample2] + 1
                        snp_dist_dict[sample2][sample1] = snp_dist_dict[sample1][sample2]
                        mlst[sample1][sample2] = mlst[sample1][sample2] + 1
                    elif s1_genotype != 'N' and s2_genotype != 'N':
                        snp_dist_rework_dict[sample1][sample2] = snp_dist_rework_dict[sample1][sample2] + 1
                        snp_dist_rework_dict[sample2][sample1] = snp_dist_rework_dict[sample1][sample2]
                        mlst_rework[sample1][sample2] = mlst[sample1][sample2] + 1
                # else:
                #	print "%s %s %s %i: N genotype excluded in SNP count" % (sample1, sample2, chrom, record.POS)

        # contig distance (unscoop to snp distance)
        for sample1, sample2 in combinations(sample_names, 2):
            if mlst[sample1][sample2] > 0:
                contig_dist_dict[sample1][sample2] = contig_dist_dict[sample1][sample2] + 1
                contig_dist_dict[sample2][sample1] = contig_dist_dict[sample1][sample2]

            if mlst_rework[sample1][sample2] > 0:
                contig_dist_rework_dict[sample1][sample2] = contig_dist_rework_dict[sample1][sample2] + 1
                contig_dist_rework_dict[sample2][sample1] = contig_dist_rework_dict[sample1][sample2]

    return snp_dist_dict, contig_dist_dict, snp_dist_rework_dict, contig_dist_rework_dict


def write_dist(dist_dic, dist_type, out_prefix):
    df = pd.DataFrame.from_records(dist_dic)
    outfile = out_prefix + '_%s_dist.tsv' % dist_type
    df.to_csv(outfile, sep='\t')
    print("{0} distance matrix written in {1}".format(dist_type, outfile))


def pre_main(args):
    vcf_file = os.path.abspath(args.vcf)
    main(vcf_file)


def main(vcf_file):
    out_prefix = os.path.splitext(vcf_file)[0]

    vcf_dic, sample_names = read_vcf(vcf_file)
    snp_dist_dict, contig_dist_dict, snp_dist_rework_dict, contig_dist_rework_dict = vcf2dist(vcf_dic, sample_names)

    write_dist(snp_dist_dict, 'SNP', out_prefix)
    write_dist(contig_dist_dict, 'Contig', out_prefix)
    write_dist(snp_dist_rework_dict, 'SNP_rework', out_prefix)
    write_dist(contig_dist_rework_dict, 'Contig_rework', out_prefix)
    print('\nvcf to matrix distances done!\n')


def version():
    return '1.0'


def run():
    parser = argparse.ArgumentParser(description='vcf2dist - Version ' + version())
    parser.add_argument('-v', '--vcf', dest="vcf", default='output_filtered.vcf', help='vcf file [output_filtered.vcf]')
    parser.add_argument('-V', '--version', action='version', version='vcf2dist-' + version(),
                        help="Prints version number")
    args = parser.parse_args()
    pre_main(args)


if __name__ == '__main__':
    run()
