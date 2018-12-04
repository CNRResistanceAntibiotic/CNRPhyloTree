#!/usr/bin/python3

import argparse
import csv
import multiprocessing
import os
import subprocess
import sys
from csv import DictReader

from skbio import DistanceMatrix
from skbio.tree import nj

from cnr_phylo_tree_src import filter_SNP_density, vcf2dist, mtx2mst


def load_config(config_file):
    config_list = []

    with open(config_file, "r") as conf:
        reader = DictReader(conf, delimiter=",")
        config_list = list(reader)

    return config_list


def run_snippy(snippy_exe, threads, out_dir, ref_genome, r1_seq_file, r2_seq_file):

    cmd = "{0} --cpus {1} --outdir {2} --reference {3} --R1 {4} --R2 {5} ".format(snippy_exe, threads, out_dir, ref_genome,
                                                                              r1_seq_file, r2_seq_file)
    log_message = " ".join(cmd)
    print(cmd)
    p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = p.communicate()
    log_message = log_message + '\n' + out.decode("utf-8") + '\n' + err.decode("utf-8")

    print(log_message)

    return log_message


def run_snippy_core_custom(snippy_exe, ref_genome, prefix, snippy_folder):
    cmd = '{0} --ref {1} --prefix {2} {3} '.format(snippy_exe, ref_genome, prefix, snippy_folder)
    print(cmd)
    os.system(cmd)


def get_snippy_dir(geno_ref_dir, result_dir, config_list):
    snippy_dir_dict = {}

    for row in config_list:
        genome_name = row["genomes"].split(".")[0]
        out_dir_root = os.path.join(result_dir, genome_name)
        list_file = os.listdir(out_dir_root)

        for file in list_file:
            file_path = os.path.join(out_dir_root, file)
            if row["strains"] in file and os.path.isfile(file_path):
                out_dir = file_path
        
        if not os.path.exists(os.path.dirname(out_dir)):
            os.mkdir(os.path.dirname(out_dir))
        ref_genome = os.path.join(geno_ref_dir, row["genomes"])

        if ref_genome not in snippy_dir_dict:
            snippy_dir_dict[ref_genome] = [out_dir]
        else:
            out_dir_list = snippy_dir_dict[ref_genome]

            out_dir_list.append(out_dir)

            # update dict
            snippy_dir_dict[ref_genome] = out_dir_list

    return snippy_dir_dict


def manage_snippy(read_dir, geno_ref_dir, result_dir, config_list):
    print("\n")

    snippy_exe = "/usr/local/snippy/bin/snippy"
    threads = 8

    strain_list = os.listdir(read_dir)

    snippy_dir_dict = {}
    jobs = []

    for row in config_list:
        genome_name = row["genomes"].split(".")[0]
        out_dir = os.path.join(result_dir, genome_name, row['strains'])

        r1_seq_file = ""
        r2_seq_file = ""

        ref_genome = os.path.join(geno_ref_dir, row["genomes"])
        if not os.path.exists(ref_genome):
            print("your reference genome file {0} don't exists or we haven't the permission\n".format(ref_genome))
            print(usage)
            sys.exit()

        for strain in strain_list:
            if row['strains']+"_" in strain:
                if "R1" in strain:
                    r1_seq_file = os.path.join(read_dir, strain)

                    if not os.path.exists(r1_seq_file):
                        print("the fastq file {0} don't exists or we haven't the permission\n".format(
                            r1_seq_file))
                        print(usage)
                        sys.exit()

                if "R2" in strain:
                    r2_seq_file = os.path.join(read_dir, strain)

                    if not os.path.exists(r2_seq_file):
                        print("the fastq file {0} don't exists or we haven't the permission\n".format(
                            r2_seq_file))
                        print(usage)
                        sys.exit()

        if not r1_seq_file:
            print("The sequence file {0} is not found".format(r1_seq_file))
            print("The program continue but without the strain  {0}".format(row['strains']))
            continue

        if not r2_seq_file:
            print("The sequence file {0} is not found".format(r1_seq_file))
            print("The program continue but without the strain  {0}".format(row['strains']))
            continue

        if ref_genome not in snippy_dir_dict:
            snippy_dir_dict[ref_genome] = [out_dir]
        else:
            out_dir_list = snippy_dir_dict[ref_genome]

            out_dir_list.append(out_dir)

            # update dict
            snippy_dir_dict[ref_genome] = out_dir_list

        if not os.path.exists(out_dir):

            p = multiprocessing.Process(
                target=run_snippy,
                args=(snippy_exe, threads, out_dir, ref_genome, r1_seq_file, r2_seq_file),
                name='snippy {0}'.format(row['strains'])
            )
            jobs.append(p)

        else:
            print("The snippy folder {0} already exist !".format(out_dir))

    job_list = []
    threshold = 3
    count = 1

    # launch and wait multiple jobs
    for job in jobs:
        job_list.append(job)
        if count == threshold:
            # launch jobs
            for job_ready in job_list:
                # remove element
                jobs.remove(job_ready)
                # start
                job_ready.start()
            # wait jobs
            for job_ready in job_list:
                job_ready.join()
                job_list = []
                count = 0
        count += 1

    # launch rest of jobs
    for job in jobs:
        job.start()

    # wait rest of jobs
    for job in jobs:
        job.join()

    return snippy_dir_dict


def manage_snippy_core(snippy_dir_dict):
    snippy_exe = "/usr/local/snippy/bin/snippy-core"
    name_dir = "core_genome"
    vcf_list = []
    print("\n")

    for genome_ref, snippy_dir_list in snippy_dir_dict.items():
        snippy_dirs = " ".join(snippy_dir_list)

        # case with only one sample that is not respectable !
        if len(snippy_dir_list) == 1:
            continue

        core_genome_path = os.path.join(os.path.dirname(snippy_dir_list[0]), name_dir)

        if not os.path.exists(core_genome_path):
            os.mkdir(core_genome_path)
            prefix_snippy_core = os.path.join(core_genome_path, name_dir)
            run_snippy_core_custom(snippy_exe, genome_ref, prefix_snippy_core, snippy_dirs)
        else:
            print("The core genome folder produces by snippy-core already exist : {0}".format(core_genome_path))

    for file in os.listdir(core_genome_path):
        if "{0}.vcf".format(name_dir) in file:
            vcf_list.append(os.path.join(core_genome_path, file))

    return vcf_list


def manage_filter_snp(vcf_list):
    min_dist = 20
    out_prefix = "filter"
    filter_keep_vcf_list = []

    for vcf_core_file in vcf_list:

        if not os.path.exists(os.path.join(os.path.dirname(vcf_core_file),
                                           out_prefix + "_" +
                                           os.path.basename(vcf_core_file).split(".vcf")[0] +
                                           "_density_filtered_keep.vcf")):

            filter_SNP_density.main(min_dist, vcf_core_file, out_prefix)

        else:
            print("the filtration is already done for {0}".format(vcf_core_file))

        for file in os.listdir(os.path.dirname(vcf_core_file)):

            if "_filtered_keep.vcf" in file:
                filter_keep_vcf_list.append(os.path.join(os.path.dirname(vcf_core_file), file))

    return filter_keep_vcf_list


def manage_r_matrix(filter_keep_vcf_list):
    r_matrix_list = []

    jump = False

    # if the matrix file are already present
    for filter_keep_vcf in filter_keep_vcf_list:

        for file in os.listdir(os.path.dirname(filter_keep_vcf)):

            if "_filtered_keep_SNP_dist.tsv" in file:
                r_matrix_list.append(os.path.join(os.path.dirname(filter_keep_vcf), file))
                jump = True
                print("the matrix step is already done for {0}".format(filter_keep_vcf))

    if not jump:

        for filter_keep_vcf in filter_keep_vcf_list:
            vcf2dist.main(filter_keep_vcf)

            for file in os.listdir(os.path.dirname(filter_keep_vcf)):

                if "_filtered_keep_SNP_dist.tsv" in file:
                    r_matrix_list.append(os.path.join(os.path.dirname(filter_keep_vcf), file))

    return r_matrix_list


def manage_make_tree(filter_keep_vcf_list):

    for filter_keep_vcf in filter_keep_vcf_list:

        filename_newick = "newick_tree.nwk"
        file_newick_path = os.path.join(os.path.dirname(filter_keep_vcf), filename_newick)

        if not os.path.exists(file_newick_path):

            with open(filter_keep_vcf, "r") as conf:
                reader = DictReader(conf, delimiter="\t")
                headers = reader.fieldnames

                # remove empty col
                headers = list(filter(None, headers))

                matrix = []
                for row in reader:
                    line = []
                    for col in headers:

                        value = row[col]
                        line.append(value)
                    matrix.append(line)

                dm = DistanceMatrix(matrix, headers)

                # print(dm)

                tree = nj(dm)

                # print(tree.ascii_art())

                newick_str = nj(dm, result_constructor=str)

                with open(file_newick_path, 'w') as out:
                    out.write("{0}".format(newick_str))

                print("the newick generation step is done for {0}".format(filter_keep_vcf))

        else:
            print("the newick generation step is already done for {0}".format(filter_keep_vcf))


def manage_snp_network(r_matrix_list, config_file, filter_keep_vcf_list):
    for mtx_file in r_matrix_list:
        graph_name = os.path.join(os.path.dirname(mtx_file), "SNP_network.html")

        with open(filter_keep_vcf_list[0]) as f:
            count_snp_keep = sum(1 for line in f)-6

        mtx2mst.main(mtx_file, graph_name, config_file, count_snp_keep)


def pre_main(args):
    read_dir = ""
    geno_ref_dir = ""
    result_dir = ""
    config_file = ""
    jump_snippy_detection = False

    if args.repRead:
        read_dir = os.path.abspath(args.repRead)

    if args.genomeRef:
        geno_ref_dir = os.path.abspath(args.genomeRef)

    if args.repResult:
        result_dir = os.path.abspath(args.repResult)

    if args.config:
        config_file = os.path.abspath(args.config)

    if args.Already:
        jump_snippy_detection = args.Already

    if not jump_snippy_detection:
        if not os.path.exists(read_dir):  # if the folder of reads exist, continue
            print("your path of directory reads don't exists or we haven't the permission\n")
            print(usage)
            sys.exit()

    if not os.path.exists(geno_ref_dir):  # if the folder of genome exist, continue
        print("your reference genome directory don't exists or we haven't the permission\n")
        print(usage)
        sys.exit()

    elif not os.path.isdir(geno_ref_dir):
        print("""\nYou have enter a file and not the folder which contain genome(s),
                with the option configfile it isn't right\n""")
        sys.exit()

    if not os.path.exists(config_file):  # if the folder of genome exist, continue
        print("your config file don't exists or we haven't the permission\n")
        print(usage)
        sys.exit()

    if not os.path.exists(result_dir):
        print("Result directory was created at '{0}'".format(result_dir))
        os.makedirs(result_dir)

    main(read_dir, geno_ref_dir, result_dir, config_file, jump_snippy_detection)


def main(read_dir, geno_ref_dir, result_dir, config_file, jump_snippy_detection=False):
    """
    This make a configuration file with argument for launch of snakemake
    :param read_dir:
    :param geno_ref_dir:
    :param result_dir:
    :param config_file:
    """

    config_list = load_config(config_file)
    snippy_dir_dict = {}


    if not jump_snippy_detection:
        # snippy
        print("\nStart Snippy")
        snippy_dir_dict = manage_snippy(read_dir, geno_ref_dir, result_dir, config_list)
        print("End Snippy")
        print("*********************************************")

    else:
        print("Skip Snippy Detection")
        snippy_dir_dict = get_snippy_dir(geno_ref_dir, result_dir, config_list)
        print("*********************************************")

    # snippy_core
    print("\nStart Snippy-core")
    vcf_list = manage_snippy_core(snippy_dir_dict)
    print("End Snippy-core")
    print("*********************************************")

    # filter SNP
    print("\nStart filtering SNP")
    filter_keep_vcf_list = manage_filter_snp(vcf_list)
    print("End filtering SNP")
    print("*********************************************")

    # Rmatrix
    print("\nStart distance matrix")
    r_matrix_list = manage_r_matrix(filter_keep_vcf_list)
    print("End distance matrix")
    print("*********************************************")


    # Make tree
    print("\nStart make tree")
    manage_make_tree(r_matrix_list)
    print("End make tree")
    print("*********************************************")

    # networkX
    print("\nStart networkX")
    manage_snp_network(r_matrix_list, config_file, filter_keep_vcf_list)
    print("End networkX")


def run():
    """
    Take the argument of the command and give a variable for this
    :return: args
    """

    global usage

    usage = "CNR Phylo tree is a snakemake based programm"

    parser = argparse.ArgumentParser(
        description="This program make a phylogenetic tree with NGS sequence reads and one or more reference genome")
    parser.add_argument('-i', '--repRead', dest="repRead", default='',
                        help='Enter the path of the directory which contain your fasta/fastq/fas/fa(.gz) reads')
    parser.add_argument('-o', '--repResult', dest="repResult", default='',
                        help="""Enter the path of your directory where you want your result (without other result or
                        risk of nothing that will be done)""")
    parser.add_argument('-g', '--repGenomeRef', dest="genomeRef", default='',
                        help="""Enter the path of directory with reference genome
                        (if no configfile enter the path to your reference genome)""")
    parser.add_argument('-c', '--config', dest="config", default='', help="Enter the path of your csv file which "
                                                                          "config the association of genome with "
                                                                          "strain (help of csv2json.py for see how "
                                                                          "make this csv file)")
    parser.add_argument('-p', '--SnippyAlreadyDone', dest="Already", action='store_true', default=False,
                        help="Indicate if all snippy detection are already done (default: False)")

    args = parser.parse_args()
    pre_main(args)


def version():
    return "1.0"


if __name__ == '__main__':
    run()
