#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
    This script will manage the vcf filtering
"""

import argparse
import os
from collections import OrderedDict

import vcf
import datetime


def read_vcf(vcf_file, not_snp_dict):
    print(f"\nRead vcf file {vcf_file}")
    vcf_dic = OrderedDict()
    vcf_reader = vcf.Reader(open(vcf_file, 'r'))
    sample_names = vcf_reader.samples
    for record in vcf_reader:
        if record.CHROM not in vcf_dic.keys():
            vcf_dic[record.CHROM] = [record]
        else:
            vcf_dic[record.CHROM].append(record)
        if f"{record.CHROM}|{record.POS}" in not_snp_dict:
            continue
        elif "snp" not in record.INFO['TYPE']:
            not_snp_dict[f"{record.CHROM}|{record.POS}"] = ""
    print(f"Number of loci   : {len(vcf_dic)}")
    print(f"Number of samples: {len(sample_names)}")
    return vcf_dic, sample_names, not_snp_dict


def pre_filter_vcf(vcf_file):
    print(f"\nRead vcf file {vcf_file}")
    vcf_reader = vcf.Reader(open(vcf_file, 'r'))
    sample_names = vcf_reader.samples

    vcf_keep_list = []
    vcf_unkeep_list = []
    count = 0

    # read line of the cvf object
    for var in vcf_reader:
        count += 1
        var_list = [var.CHROM, var.POS, '.', var.REF, var.ALT, ".", "PASS", var.INFO, var.FORMAT]
        snp_name_hash = {}

        for sample in sample_names:
            snp_name_hash[sample] = var.genotype(sample)['GT']
        var_list.append(snp_name_hash)

        res_set = set(list(var_list[9].values()))

        if len(res_set) == 1:
            vcf_unkeep_list.append(var_list)
        elif "N" in var_list[4] and len(res_set) == 2:
            vcf_unkeep_list.append(var_list)
        else:
            vcf_keep_list.append(var_list)

    return vcf_keep_list, vcf_unkeep_list, sample_names, count


def write_vcf(vcf_list, sample_names, out_vcf_file):
    """
    This function write a vcf file by a list of snps in vcf format
    :param vcf_list: a list of snps
    :param sample_names: the name of sample
    :param out_vcf_file: the path of the output file
    :return: nothing
    """
    now = datetime.datetime.now()
    col_names = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT']
    header = "##fileformat=VCFv4.2\n"
    header = header + f"##fileDate={now.year}{now.month}{now.day}\n"
    header = header + "##source=filter_SNP_densityV0.1\n"
    header = header + "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
    header = header + "##INFO=<ID=TYPE,Number=A,Type=String,Description=\"Allele type: snp ins del\">\n"
    header = "{0}{1}".format(header, "#{0}\t{1}\n".format('\t'.join(col_names), '\t'.join(sample_names)))
    with open(out_vcf_file, 'w') as vcf_output_file:
        vcf_output_file.write(header)
        for element in vcf_list:
            line_to_write = ""
            count = 0
            for x in element:
                count += 1
                # final dict
                if count == 10:
                    for sample in sample_names:
                        line_to_write = line_to_write + str(x.get(sample)) + "\t"
                elif type(x) == list:
                    line_to_write = line_to_write + str(x).replace("[", "").replace("]", "").replace("'", "").replace(
                        " ",
                        "") + "\t"
                elif type(x) == dict:
                    # get first key
                    key = next(iter(x))
                    line_to_write = line_to_write + key + "=" + str(x.get(key)).replace("[", "").replace("]",
                                                                                                         "").replace(
                        "'", "") + "\t"
                else:
                    line_to_write = line_to_write + str(x).replace("'", "").replace(" ", "") + "\t"
            line_to_write += "\n"
            vcf_output_file.write(line_to_write)


def pre_main(args):
    vcf_file = os.path.abspath(args.vcf)
    out_vcf_file = args.out_vcf_file
    out_vcf_unkeep_file = args.out_vcf_unkeep_file
    main(vcf_file, out_vcf_file, out_vcf_unkeep_file)


def main(vcf_file, out_vcf_file, out_vcf_unkeep_file):

    print("\nDATE: ", datetime.datetime.now())

    vcf_list, vcf_unkeep_list, sample_names, count = pre_filter_vcf(vcf_file)
    print("\nDATE: ", datetime.datetime.now())
    vcf_list = sorted(vcf_list, key=lambda element: (element[0], element[1]))
    vcf_unkeep_list = sorted(vcf_unkeep_list, key=lambda element: (element[0], element[1]))
    print("\n " + str(count) + " vcf entries are read.")
    print(f"\n             {len(vcf_list)} filtered vcf entries are saved. {round((len(vcf_list) * 100) / count)}%")
    print(f"\n             {len(vcf_unkeep_list)} filtered vcf entries are saved in unkeep file. {round((len(vcf_unkeep_list) * 100) / count)}%")
    if len(vcf_list) + len(vcf_unkeep_list) == count:
        print("\n SUCCESS to pre-filter all SNPs")
    else:
        print("\n FAIL to pre-filter all SNPs")

    write_vcf(vcf_list, sample_names, out_vcf_file)
    write_vcf(vcf_unkeep_list, sample_names, out_vcf_unkeep_file)

    print("\nDATE: ", datetime.datetime.now())
    print('\nVCF pre-filtration done!\n')


def version():
    return '1.0'


def run():
    parser = argparse.ArgumentParser(description='VCFPreFilterSNPDensity- Version ' + version())
    parser.add_argument('-v', '--vcf', dest="vcf", default='test.vcf', help='vcf file [snps.recode.vcf]')
    parser.add_argument('-o', '--out', dest="out_vcf_file", default='', help='output kepp')
    parser.add_argument('-ou', '--outunkeep', dest="out_vcf_unkeep_file", default='', help='output unkeep')
    parser.add_argument('-V', '--version', action='version', version='VCFPreFilterSNPDensity-' + version(),
                        help="Prints version number")
    args = parser.parse_args()
    pre_main(args)


if __name__ == '__main__':
    run()
