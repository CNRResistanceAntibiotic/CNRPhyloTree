#!/usr/bin/python3
# coding=utf-8

import argparse
import os
import vcf
import datetime


def read_vcf(vcf_file, threshold=10):
    print("\nRead vcf file {}".format(vcf_file))

    vcf_reader = vcf.Reader(open(vcf_file, 'r'))
    sample_names = vcf_reader.samples

    old_var_chrom = ""
    vcf_keep_list = []
    vcf_unkeep_list = []
    count = 0
    snp_hash = {}

    # read line of the cvf object
    for var in vcf_reader:

        count += 1

        var_list = [var.CHROM, var.POS, '.', var.REF, var.ALT, var.QUAL, '.', var.INFO, var.FORMAT]

        snp_name_hash = {}
        for sample in sample_names:
            snp_name_hash[sample] = var.genotype(sample)['GT']

        var_list.append(snp_name_hash)

        #  case of first snp of the first CHROM
        if not snp_hash:
            snp_hash[var.POS] = var_list
            old_var_chrom = var.CHROM
            continue

        else:

            if str(old_var_chrom) == str(var.CHROM):
                snp_hash[var.POS] = var_list;
                continue

            else:

                snp_hash = update_vcf_hash(snp_hash, sample_names, threshold)

                # class snp to keep or unkeep
                for pos in snp_hash:
                    keep_list, unkeep_list = class_vcf(snp_hash.get(pos))
                    if keep_list:
                        vcf_keep_list.append(keep_list)
                    if unkeep_list:
                        vcf_unkeep_list.append(unkeep_list)

                # reset hash for next chrom
                snp_hash.clear()

                # first snp of the actual new chrom
                snp_hash[var.POS] = var_list

                # new chrom ref
                old_var_chrom = var.CHROM

    #######################################################################
    # last chrom
    snp_hash = update_vcf_hash(snp_hash, sample_names, threshold)

    # class snp to keep or unkeep
    for pos in snp_hash:
        keep_list, unkeep_list = class_vcf(snp_hash.get(pos))
        if keep_list:
            vcf_keep_list.append(keep_list)
        if unkeep_list:
            vcf_unkeep_list.append(unkeep_list)


    return vcf_keep_list, vcf_unkeep_list, sample_names, count;


def update_vcf_hash(snp_hash, sample_names, threshold):
    """
    This function update the vcf hash
    :param snp_hash: a vcf of snps hash
    :param sample_names: the list of sample names
    :param threshold: the threshold of min distance btw 2 snp
    :return: a hash update of snp for 1 chrom
    """
    for pos_1 in snp_hash:
        for pos_2 in snp_hash:

            # dumb case
            if pos_1 is not pos_2:
                # make the threshold test
                if abs(pos_1 - pos_2) <= threshold:

                    var_list_1 = snp_hash.get(pos_1)
                    var_list_2 = snp_hash.get(pos_2)

                    for sample in sample_names:

                        # compare base of 2 snp for the same sample
                        if var_list_1[9][sample] is not "0":
                            if var_list_2[9][sample] is not "0":

                                # add "N" to the ALT
                                if "N" not in var_list_1[4]:
                                    var_list_1[4].append("N")

                                # add "N" to the ALT
                                if "N" not in var_list_2[4]:
                                    var_list_2[4].append("N")

                                # update number ALT
                                var_list_1[9][sample] = len(var_list_1[4])
                                var_list_2[9][sample] = len(var_list_2[4])

                                # update hash
                                snp_hash[pos_1] = var_list_1
                                snp_hash[pos_2] = var_list_2
    return snp_hash


def class_vcf(var_list):
    """
    This function classify the snp in two list : keep and unkeep
    :param var_list: a snp in list format.
    :return: two lists : keep and unkeep
    """

    vcf_keep_list = []
    vcf_unkeep_list = []

    # get list of snp value for all samples
    res_set = set(list(var_list[9].values()))

    # case with only ref and undetermined allele is found on vcf entry -> unkeep
    if "N" in var_list[4] and "0" in res_set and len(var_list[4]) in res_set and len(res_set) == 2:
        vcf_unkeep_list = var_list

    # case with only the same allele is found on vcf entry -> unkeep
    elif len(res_set) == 1:
        vcf_unkeep_list = var_list

    # case with only the undetermined allele is found on vcf entry -> unkeep
    elif "N" in var_list[4] and len(var_list[4]) in res_set and len(res_set) == 1:
        vcf_unkeep_list = var_list

    # case with only the ref allele is found on vcf entry-> unkeep
    elif "0" in res_set and len(res_set) == 1:
        vcf_unkeep_list = var_list

    # all others cases -> keep
    else:
        vcf_keep_list = var_list

    return vcf_keep_list, vcf_unkeep_list


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

    header = "##fileformat=VCFv4.1\n"
    header = header + "##fileDate={0}{1}{2}\n".format(now.year, now.month, now.day)
    header = header + "##source=filter_SNP_densityV0.1\n"
    header = header + "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
    header = header + "##INFO=<ID=TYPE,Number=A,Type=String,Description=\"The type of allele.\">\n"
    header = header + "#{0}\t{1}\n".format('\t'.join(col_names), '\t'.join(sample_names))
    with open(out_vcf_file, 'w') as vcf_output_file:
        vcf_output_file.write(header)
        for element in vcf_list:

            line_to_write = ""
            count = 0;
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
    min_dist = int(args.mindist)
    vcf_file = os.path.abspath(args.vcf)
    out_prefix = args.outPrefix
    main(min_dist, vcf_file, out_prefix)


def main(min_dist, vcf_file, out_prefix):

    print("\nSNP density threshold: 2 SNPs for {0} base(s)".format(min_dist))

    dir_name_path = os.path.dirname(vcf_file)
    name_file_vcf_input = os.path.basename(os.path.splitext(vcf_file)[0])

    if out_prefix:
        name_output_file = out_prefix + "_" + name_file_vcf_input
        out_vcf_file = os.path.join(dir_name_path, name_output_file + "_density_filtered_keep.vcf")
        out_vcf_unkeep_file = os.path.join(dir_name_path, name_output_file + "_density_filtered_unkeep.vcf")
    else:
        out_vcf_file = os.path.join(dir_name_path, name_file_vcf_input + "_density_filtered_keep.vcf")
        out_vcf_unkeep_file = os.path.join(dir_name_path, name_file_vcf_input + "_density_filtered_unkeep.vcf")

    print("\nDATE: ", datetime.datetime.now())

    vcf_list, vcf_unkeep_list, sample_names, count = read_vcf(vcf_file, min_dist)

    print("\nDATE: ", datetime.datetime.now())

    vcf_list = sorted(vcf_list, key=lambda element: (element[0], element[1]))
    vcf_unkeep_list = sorted(vcf_unkeep_list, key=lambda element: (element[0], element[1]))

    print("\n " + str(count) + " vcf entries are read.")
    print("\n             {0} filtered vcf entries are saved. {1}%"
          .format(len(vcf_list), round((len(vcf_list)*100)/count) ,2))
    print("\n             {0} filtered vcf entries are saved in unkeep file. {1}%"
          .format(len(vcf_unkeep_list), round((len(vcf_unkeep_list)*100)/count),2))
    if len(vcf_list) + len(vcf_unkeep_list) == count:
        print("\n SUCCESS to class all SNPs")
    else:
        print("\n FAIL to class all SNPs")

    write_vcf(vcf_list, sample_names, out_vcf_file)
    write_vcf(vcf_unkeep_list, sample_names, out_vcf_unkeep_file)

    print("\nDATE: ", datetime.datetime.now())

    print('\nvcf filtration done!\n')


def version():
    return '1.0'


def run():
    parser = argparse.ArgumentParser(description='VCFfilterSNPdensity- Version ' + version())
    parser.add_argument('-v', '--vcf', dest="vcf", default='test.vcf', help='vcf file [snps.recode.vcf]')
    parser.add_argument('-min', '--mindist', dest="mindist", default='20',
                        help='Select minimum distance between SNPs [10]')
    parser.add_argument('-o', '--out', dest="outPrefix", default='', help='ouput prefix')
    parser.add_argument('-V', '--version', action='version', version='VCFfilterSNPdensity-' + version(),
                        help="Prints version number")
    args = parser.parse_args()
    pre_main(args)


if __name__ == '__main__':
    run()
