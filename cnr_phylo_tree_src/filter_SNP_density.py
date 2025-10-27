#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
VCF filtering script (compatible with vcfpy)
"""

import argparse
import os
import datetime
import vcfpy


def read_with_min_dist_vcf(vcf_file, threshold, type_matrix):
    print(f"\nRead vcf file {vcf_file}")
    vcf_reader = vcfpy.Reader.from_path(vcf_file)
    sample_names = vcf_reader.header.samples.names
    old_var_chrom = ""
    vcf_keep_list, vcf_unkeep_list = [], []
    count = 0
    snp_hash = {}

    for var in vcf_reader:
        count += 1
        var_list = [
            var.CHROM,
            var.POS,
            ".",
            var.REF,
            [alt.value for alt in var.ALT],
            var.QUAL,
            ".",
            var.INFO,
            ":".join(var.FORMAT)
        ]

        # dictionnaire des génotypes par échantillon
        snp_name_hash = {}
        for sample in sample_names:
            call = var.call_for_sample.get(sample)
            snp_name_hash[sample] = call.data.get("GT") if call else None
        var_list.append(snp_name_hash)

        if not snp_hash:
            snp_hash[var.POS] = var_list
            old_var_chrom = var.CHROM
            continue
        else:
            if str(old_var_chrom) == str(var.CHROM):
                snp_hash[var.POS] = var_list
                continue
            else:
                snp_hash = update_vcf_hash(snp_hash, sample_names, threshold)
                for pos in snp_hash:
                    keep_list, unkeep_list = class_vcf(snp_hash.get(pos), type_matrix)
                    if keep_list:
                        vcf_keep_list.append(keep_list)
                    if unkeep_list:
                        vcf_unkeep_list.append(unkeep_list)
                snp_hash.clear()
                snp_hash[var.POS] = var_list
                old_var_chrom = var.CHROM

    # dernier chromosome
    print("last chrom")
    print(snp_hash, sample_names, threshold)
    snp_hash = update_vcf_hash(snp_hash, sample_names, threshold)
    print(snp_hash)
    print("\nClass SNPs in Keep or UnKeep\n")
    for pos in snp_hash:
        keep_list, unkeep_list = class_vcf(snp_hash.get(pos), type_matrix)
        if keep_list:
            vcf_keep_list.append(keep_list)
        if unkeep_list:
            vcf_unkeep_list.append(unkeep_list)

    return vcf_keep_list, vcf_unkeep_list, sample_names, count


def update_vcf_hash(snp_hash, sample_names, threshold):
    for pos_1 in snp_hash:
        for pos_2 in snp_hash:
            if pos_1 != pos_2:
                if abs(int(pos_1) - int(pos_2)) <= int(threshold):
                    var_list_1 = snp_hash.get(pos_1)
                    var_list_2 = snp_hash.get(pos_2)
                    for sample in sample_names:
                        if var_list_1[9][sample] != "0" and var_list_2[9][sample] != "0":
                            if "X" not in var_list_1[4]:
                                var_list_1[4].append("X")
                            if "X" not in var_list_2[4]:
                                var_list_2[4].append("X")
                            var_list_1[9][sample] = len(var_list_1[4])
                            var_list_2[9][sample] = len(var_list_2[4])
                            snp_hash[pos_1] = var_list_1
                            snp_hash[pos_2] = var_list_2
    return snp_hash


def class_vcf(var_list, type_matrix):
    vcf_keep_list, vcf_unkeep_list = [], []
    res_set = set(list(var_list[9].values()))

    if "X" in var_list[4]:
        vcf_unkeep_list = var_list
    elif "N" in var_list[4]:
        if type_matrix == "absolu":
            vcf_unkeep_list = var_list
        elif type_matrix in ["relatif", "relatif-norm"]:
            vcf_keep_list = var_list
    elif len(res_set) == 1:
        vcf_unkeep_list = var_list
    else:
        vcf_keep_list = var_list
    return vcf_keep_list, vcf_unkeep_list


def write_vcf(vcf_list, sample_names, out_vcf_file):
    now = datetime.datetime.now()
    col_names = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT']
    header = f"##fileformat=VCFv4.2\n##fileDate={now.year}{now.month}{now.day}\n"
    header += "##source=filter_SNP_density_vcfpy\n"
    header += "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
    header += "##INFO=<ID=TYPE,Number=A,Type=String,Description=\"The type of allele.\">\n"
    header += "#{0}\t{1}\n".format('\t'.join(col_names), '\t'.join(sample_names))

    with open(out_vcf_file, 'w') as vcf_output_file:
        vcf_output_file.write(header)
        for element in vcf_list:
            line_to_write = ""
            count = 0
            for x in element:
                count += 1
                if count == 10:
                    for sample in sample_names:
                        line_to_write += str(x.get(sample)) + "\t"
                elif isinstance(x, list):
                    x = ','.join(map(str, x))
                    line_to_write += x + "\t"
                elif isinstance(x, dict):
                    if len(x) == 1:
                        k = next(iter(x))
                        line_to_write += f"{k}={x[k]}\t"
                    else:
                        line_to_write += str(x) + "\t"
                else:
                    line_to_write += str(x) + "\t"
            vcf_output_file.write(line_to_write.strip() + "\n")


def pre_main(args):
    min_dist = int(args.mindist)
    vcf_file = os.path.abspath(args.vcf)
    out_prefix = args.outPrefix
    type_matrix = args.type_matrice
    main(min_dist, vcf_file, out_prefix, type_matrix)


def main(min_dist, vcf_file, out_prefix, type_matrix):
    print(f"\nSNP density threshold: 2 SNPs for {min_dist} base(s)")
    dir_name_path = os.path.dirname(vcf_file)
    name_file_vcf_input = os.path.basename(os.path.splitext(vcf_file)[0])

    if out_prefix:
        name_output_file = out_prefix + "_" + name_file_vcf_input
        out_vcf_file = os.path.join(dir_name_path, name_output_file + "_density_filtered_keep.vcf")
        out_vcf_unkeep_file = os.path.join(dir_name_path, name_output_file + "_density_filtered_unkeep.vcf")
    else:
        out_vcf_file = os.path.join(dir_name_path, name_file_vcf_input + "_density_filtered_keep.vcf")
        out_vcf_unkeep_file = os.path.join(dir_name_path, name_file_vcf_input + "_density_filtered_unkeep.vcf")

    print("\nDATE:", datetime.datetime.now())
    vcf_list, vcf_unkeep_list, sample_names, count = read_with_min_dist_vcf(vcf_file, min_dist, type_matrix)
    print("\nDATE:", datetime.datetime.now())
    vcf_list = sorted(vcf_list, key=lambda e: (e[0], e[1]))
    vcf_unkeep_list = sorted(vcf_unkeep_list, key=lambda e: (e[0], e[1]))
    print(f"\n{count} VCF entries read.")

    if count == 0:
        print("No SNP kept for further tasks")
        return False
    else:
        print(f"\n{len(vcf_list)} SNPs kept ({round(len(vcf_list)*100/count)}%)")
        print(f"{len(vcf_unkeep_list)} SNPs unkept ({round(len(vcf_unkeep_list)*100/count)}%)")
        if len(vcf_list) + len(vcf_unkeep_list) == count:
            print("\nSUCCESS: all SNPs classified")
        else:
            print("\nWARNING: not all SNPs classified")
        write_vcf(vcf_list, sample_names, out_vcf_file)
        write_vcf(vcf_unkeep_list, sample_names, out_vcf_unkeep_file)
    print("\nDATE:", datetime.datetime.now())
    print('\nVCF filtration done!\n')
    return True


def version():
    return '2.0-vcfpy'


def run():
    parser = argparse.ArgumentParser(description='VCFFilterSNPDensity - Version ' + version())
    parser.add_argument('-v', '--vcf', dest="vcf", default='test.vcf', help='vcf file [snps.recode.vcf]')
    parser.add_argument('-min', '--mindist', dest="mindist", default='20',
                        help='Select minimum distance between SNPs [10]')
    parser.add_argument('-o', '--out', dest="outPrefix", default='', help='output prefix')
    parser.add_argument('-tm', '--type_matrice', dest="type_matrice", default="absolu",
                        help="Matrix type [absolu, relatif, relatif-norm] (default: absolu)")
    parser.add_argument('-V', '--version', action='version', version='VCFFilterSNPDensity-' + version())
    args = parser.parse_args()
    pre_main(args)


if __name__ == '__main__':
    run()
