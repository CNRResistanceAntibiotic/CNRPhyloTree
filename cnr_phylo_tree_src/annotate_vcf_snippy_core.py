#!/usr/bin/python3
# -*- coding: utf-8 -*-
import argparse
import csv
import os
from collections import OrderedDict

import vcf


def read_vcf(vcf_file):
    print("\nRead vcf file {}".format(vcf_file))
    vcf_dic = OrderedDict()
    vcf_reader = vcf.Reader(open(vcf_file, 'r'))
    sample_names = vcf_reader.samples
    counter = 0
    for record in vcf_reader:
        counter += 1
        if record.CHROM not in vcf_dic.keys():
            vcf_dic[record.CHROM] = [record]
        else:
            vcf_dic[record.CHROM].append(record)
    print("Number of loci   : {}".format(len(vcf_dic)))
    print("Number of samples: {}".format(len(sample_names)))
    return vcf_dic, sample_names, counter


def annotate_vcf(vcf_dic, sample_names, vcf_folder_list, annotate_vcf_snippy_core):
    out_dict_list = []
    counter = 0
    for chrom in vcf_dic.keys():
        counter += 1
        records_filt = vcf_dic[chrom]
        for vcf_folder in vcf_folder_list:
            vcf_file = os.path.join(vcf_folder, "snps.vcf")
            vcf_strain_dict, sample_names_s, count = read_vcf(vcf_file)
            strain_name = sample_names_s[0]
            print(chrom)
            if chrom in vcf_strain_dict:
                print(vcf_strain_dict[chrom])
                records = vcf_strain_dict[chrom]
                for record_filt in records_filt:
                    for record in records:
                        if record_filt.CHROM == record.CHROM and record_filt.POS == record.POS:
                            records.remove(record)
                            out_dict = {}
                            existed = False
                            strain_name_split_list = strain_name.split("_")
                            strain_name_clean = strain_name_split_list[2]
                            for out_dict_l in out_dict_list:
                                if out_dict_l:
                                    if out_dict_l["CHROM"] == record_filt.CHROM and out_dict_l["POS"] == record_filt.POS:
                                        existed = True
                                        out_dict = out_dict_l
                                        continue
                            if not existed:
                                if "ANN" in record.INFO:
                                    out_dict["ANN"] = record.INFO["ANN"]
                                out_dict["CHROM"] = record_filt.CHROM
                                out_dict["POS"] = record_filt.POS
                                out_dict["REF"] = record_filt.REF
                                out_dict["ALT"] = record_filt.ALT
                            if record_filt.genotype(strain_name_clean)['GT'] == 0:
                                out_dict[strain_name_clean] = 0
                                continue
                            else:
                                out_dict[strain_name_clean] = f"DP:{record.INFO['DP']};RO={record.INFO['RO']};" \
                                                              f"RA={record.INFO['AO'][0]}({record.ALT[0]}/" \
                                                              f"{record_filt.ALT[int(record_filt.genotype(strain_name_clean)['GT']) - 1]})"
                            if out_dict:
                                if not existed:
                                    out_dict_list.append(out_dict)
                                break
                        else:
                            continue
            else:
                continue
    out_dict_list = sorted(out_dict_list, key=lambda i: (i['CHROM'], i['POS']))
    with open(annotate_vcf_snippy_core, "w") as output_file:
        tsv_writer = csv.writer(output_file, delimiter='\t')
        tsv_writer.writerow(["CHROM", "POS", "REF", "ALT", "ANN"] + sample_names)
        for out_dict in out_dict_list:
            sample_l = []
            for s in sample_names:
                if s in out_dict:
                    sample_l.append(out_dict[s])
                else:
                    sample_l.append("0")
            if "ANN" in out_dict:
                tsv_writer.writerow([out_dict["CHROM"], out_dict["POS"], out_dict["REF"], out_dict["ALT"],
                                     out_dict["ANN"]] + sample_l)
            else:
                tsv_writer.writerow(
                    [out_dict["CHROM"], out_dict["POS"], out_dict["REF"], out_dict["ALT"], ""] + sample_l)


def pre_main(args):
    vcf_file = os.path.abspath(args.vcf)
    vcf_folder = os.path.abspath(args.vcf_folder)
    annotate_vcf_snippy_core = os.path.abspath(args.output)
    main(vcf_file, vcf_folder, annotate_vcf_snippy_core)


def main(vcf_file, vcf_folder_list, annotate_vcf_snippy_core):
    # read snippy-core result vcf
    vcf_dic, sample_names, counter = read_vcf(vcf_file)
    threshold = 5000
    if counter <= threshold:
        # annotate with other vcf
        annotate_vcf(vcf_dic, sample_names, vcf_folder_list, annotate_vcf_snippy_core)
    else:
        print("The number of SNP in the core vcf is superior to {0} ({1}). No Supplementary annotation"
              " provided.".format(threshold, counter))


def version():
    return '1.0'


def run():
    parser = argparse.ArgumentParser(description='vcf2dist - Version ' + version())
    parser.add_argument('-v', '--vcf', dest="vcf", default='snps.vcf', help='vcf file [snps.vcf]')
    parser.add_argument('-vf', '--vcf_folder', dest="vcf_folder", default='vcf_output_snippy',
                        help='folder with strain snippy vcf')
    parser.add_argument('-o', '--output', dest="output", default='', help='annotated vcf output')
    parser.add_argument('-V', '--version', action='version', version='vcf2dist-' + version(),
                        help="Prints version number")
    args = parser.parse_args()
    pre_main(args)


if __name__ == '__main__':
    run()
