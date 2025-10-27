#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Annotate a core VCF file with SNPs from multiple Snippy strain VCFs.
Compatible with Python ≥ 3.10 and vcfpy ≥ 0.13.
"""

import argparse
import csv
import os
from collections import OrderedDict
import vcfpy


def read_vcf(vcf_file):
    """Read a VCF file and return an OrderedDict keyed by chromosome."""
    print(f"\nReading VCF file: {vcf_file}")
    vcf_dic = OrderedDict()
    reader = vcfpy.Reader.from_path(vcf_file)
    sample_names = reader.header.samples.names
    counter = 0

    for record in reader:
        counter += 1
        print(counter)
        print(record)

    print(f"  Number of loci   : {len(vcf_dic)}")
    print(f"  Number of samples: {len(sample_names)}")
    return vcf_dic, sample_names, counter


def annotate_vcf(vcf_dic, sample_names, vcf_folder_list, output_file):
    """Annotate a reference/core VCF with SNP data from multiple strain VCFs."""
    out_dict_list = []
    chrom_use_dict = {}

    # Read all strain VCFs
    for vcf_folder in vcf_folder_list:
        vcf_file = os.path.join(vcf_folder, "snps.vcf")
        if not os.path.exists(vcf_file):
            print(f"⚠️ Skipping {vcf_folder}: snps.vcf not found")
            continue

        vcf_strain_dict, strain_sample_names, _ = read_vcf(vcf_file)
        strain_name = strain_sample_names[0] if strain_sample_names else os.path.basename(vcf_folder)

        for chrom, records in vcf_strain_dict.items():
            chrom_use_dict.setdefault(chrom, []).append({"value": records, "strain_name": strain_name})

    # Annotate records
    for chrom, records_filt in vcf_dic.items():
        if chrom not in chrom_use_dict:
            continue

        for element in chrom_use_dict[chrom]:
            records = element["value"]
            strain_name = element["strain_name"]

            for record_filt in records_filt:
                for record in list(records):  # copy to allow safe removal
                    if record_filt.CHROM == record.CHROM and record_filt.POS == record.POS:
                        records.remove(record)
                        out_dict = next(
                            (d for d in out_dict_list if d and d["CHROM"] == record_filt.CHROM and d["POS"] == record_filt.POS),
                            None
                        )

                        if not out_dict:
                            out_dict = {
                                "CHROM": record_filt.CHROM,
                                "POS": int(record_filt.POS),
                                "REF": record_filt.REF,
                                "ALT": ",".join(str(a) for a in record_filt.ALT),
                                "ANN": record.INFO.get("ANN", "")
                            }
                            out_dict_list.append(out_dict)

                        geno = record_filt.call_for_sample.get(strain_name)
                        gt = geno.data.get("GT") if geno and geno.data else None

                        if gt in [None, "0", "0/0"]:
                            out_dict[strain_name] = 0
                        else:
                            dp = record.INFO.get("DP", "")
                            ro = record.INFO.get("RO", "")
                            ao = record.INFO.get("AO", [0])
                            out_dict[strain_name] = f"DP:{dp};RO={ro};RA={ao[0]}({record.ALT[0]})"

    # Sort output
    out_dict_list.sort(key=lambda i: (i["CHROM"], i["POS"]))

    # Write TSV output
    with open(output_file, "w", newline='') as out_f:
        writer = csv.writer(out_f, delimiter='\t')
        writer.writerow(["CHROM", "POS", "REF", "ALT", "ANN"] + sample_names)
        for out_dict in out_dict_list:
            row = [
                out_dict["CHROM"],
                out_dict["POS"],
                out_dict["REF"],
                out_dict["ALT"],
                out_dict.get("ANN", "")
            ] + [out_dict.get(s, "0") for s in sample_names]
            writer.writerow(row)

    print(f"\n✅ Annotation written to: {output_file}")


def main(vcf_file, vcf_folder_list, output_file):
    vcf_dic, sample_names, count = read_vcf(vcf_file)
    threshold = 5000
    if count <= threshold:
        annotate_vcf(vcf_dic, sample_names, vcf_folder_list, output_file)
    else:
        print(f"⚠️ The number of SNPs in the core VCF exceeds {threshold} ({count}). Annotation skipped.")


def run():
    parser = argparse.ArgumentParser(description="Annotate a core VCF with Snippy strain VCFs")
    parser.add_argument("-v", "--vcf", required=True, help="Core VCF file (e.g. snps.vcf)")
    parser.add_argument("-vf", "--vcf_folder", nargs='+', required=True,
                        help="List of folders each containing a Snippy strain VCF")
    parser.add_argument("-o", "--output", required=True, help="Output TSV file")
    args = parser.parse_args()

    vcf_file = os.path.abspath(args.vcf)
    vcf_folder_list = [os.path.abspath(f) for f in args.vcf_folder]
    output_file = os.path.abspath(args.output)

    main(vcf_file, vcf_folder_list, output_file)


if __name__ == "__main__":
    run()
