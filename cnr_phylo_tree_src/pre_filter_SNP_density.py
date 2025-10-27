#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
    This script will manage the VCF filtering
    Compatible with Python >= 3.10 using vcfpy
"""

import argparse
import os
from collections import OrderedDict
import datetime
import vcfpy


def read_vcf(vcf_file, not_snp_dict):
    print(f"\nRead VCF file {vcf_file}")
    vcf_dic = OrderedDict()
    reader = vcfpy.Reader.from_path(vcf_file)
    sample_names = reader.header.samples.names

    for record in reader:
        chrom = record.CHROM
        if chrom not in vcf_dic:
            vcf_dic[chrom] = [record]
        else:
            vcf_dic[chrom].append(record)

        # Gestion du dictionnaire des positions non SNP
        key = f"{record.CHROM}|{record.POS}"
        if key in not_snp_dict:
            continue
        elif "snp" not in record.INFO.get("TYPE", []):
            not_snp_dict[key] = ""

    print(f"Number of loci   : {len(vcf_dic)}")
    print(f"Number of samples: {len(sample_names)}")
    return vcf_dic, sample_names, not_snp_dict


def pre_filter_vcf(vcf_file):
    print(f"\nRead VCF file {vcf_file}")
    reader = vcfpy.Reader.from_path(vcf_file)
    sample_names = reader.header.samples.names

    vcf_keep_list = []
    vcf_unkeep_list = []
    count = 0

    for var in reader:
        count += 1
        # Champs du VCF
        var_list = [
            var.CHROM,
            var.POS,
            ".",
            var.REF,
            [str(alt) for alt in var.ALT],
            ".",
            "PASS",
            var.INFO,
            var.FORMAT,
        ]
        snp_name_hash = {}

        for sample in sample_names:
            gt = var.call_for_sample[sample].data.get("GT")
            snp_name_hash[sample] = gt if gt is not None else "./."

        var_list.append(snp_name_hash)
        res_set = set(var_list[9].values())

        # Conditions d'exclusion ou d'inclusion
        if len(res_set) == 1:
            vcf_unkeep_list.append(var_list)
        elif "N" in var_list[4] and len(res_set) == 2:
            vcf_unkeep_list.append(var_list)
        else:
            vcf_keep_list.append(var_list)

    return vcf_keep_list, vcf_unkeep_list, sample_names, count


def write_vcf(vcf_list, sample_names, out_vcf_file):
    """
    Write a VCF file from a list of SNP entries
    """
    now = datetime.datetime.now()
    col_names = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT']

    header = "##fileformat=VCFv4.2\n"
    header += f"##fileDate={now.year}{now.month:02d}{now.day:02d}\n"
    header += "##source=filter_SNP_densityV0.1\n"
    header += "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
    header += "##INFO=<ID=TYPE,Number=A,Type=String,Description=\"Allele type: snp ins del\">\n"
    header += "#{0}\t{1}\n".format('\t'.join(col_names), '\t'.join(sample_names))

    with open(out_vcf_file, 'w') as vcf_output_file:
        vcf_output_file.write(header)
        for element in vcf_list:
            line_to_write = ""
            count = 0
            for x in element:
                count += 1
                if count == 10:  # dictionnaire des génotypes
                    for sample in sample_names:
                        gt = x.get(sample, "./.")
                        line_to_write += gt + "\t"
                elif isinstance(x, list):
                    line_to_write += ','.join(map(str, x)) + "\t"
                elif isinstance(x, dict):
                    # Écrire sous forme clé=valeur
                    if len(x) > 0:
                        key = next(iter(x))
                        val = x[key]
                        if isinstance(val, list):
                            val = ",".join(map(str, val))
                        line_to_write += f"{key}={val}\t"
                    else:
                        line_to_write += ".\t"
                else:
                    line_to_write += str(x) + "\t"

            line_to_write = line_to_write.strip() + "\n"
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

    vcf_list = sorted(vcf_list, key=lambda e: (e[0], e[1]))
    vcf_unkeep_list = sorted(vcf_unkeep_list, key=lambda e: (e[0], e[1]))

    print(f"\n {count} VCF entries are read.")
    print(f"\n {len(vcf_list)} kept ({round((len(vcf_list) * 100) / count)}%)")
    print(f"\n {len(vcf_unkeep_list)} unkept ({round((len(vcf_unkeep_list) * 100) / count)}%)")

    if len(vcf_list) + len(vcf_unkeep_list) == count:
        print("\n SUCCESS: all SNPs pre-filtered.")
    else:
        print("\n WARNING: mismatch between total and filtered entries.")

    write_vcf(vcf_list, sample_names, out_vcf_file)
    write_vcf(vcf_unkeep_list, sample_names, out_vcf_unkeep_file)

    print("\nDATE: ", datetime.datetime.now())
    print('\nVCF pre-filtration done!\n')

    if len(vcf_list) == 0:
        print("\nNo SNP FOUND! Exit program.")
        exit(1)


def version():
    return '1.1'


def run():
    parser = argparse.ArgumentParser(description='VCFPreFilterSNPDensity - Version ' + version())
    parser.add_argument('-v', '--vcf', dest="vcf", default='test.vcf', help='VCF file [snps.recode.vcf]')
    parser.add_argument('-o', '--out', dest="out_vcf_file", default='filtered_keep.vcf', help='output keep file')
    parser.add_argument('-ou', '--outunkeep', dest="out_vcf_unkeep_file", default='filtered_unkeep.vcf',
                        help='output unkeep file')
    parser.add_argument('-V', '--version', action='version', version='VCFPreFilterSNPDensity-' + version(),
                        help="Prints version number")
    args = parser.parse_args()
    pre_main(args)


if __name__ == '__main__':
    run()
