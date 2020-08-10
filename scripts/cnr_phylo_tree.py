#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
    This script will manage the phylo tree construction
"""

import argparse
import multiprocessing
import os
import subprocess
import sys
from csv import DictReader
from xml.dom import minidom

from Bio import Phylo
# from Bio.Phylo.TreeConstruction import DistanceTreeConstructor, DistanceMatrix
from skbio import DistanceMatrix
from skbio.tree import nj

from cnr_phylo_tree_src import filter_SNP_density, vcf2dist, mtx2mst, annotate_vcf_snippy_core
import xml.etree.ElementTree as ET


def load_config(config_file):
    """
    This function load config file
    :param config_file:the config file path
    :return: the config data in list format
    """

    with open(config_file, "r") as conf:
        reader = DictReader(conf, delimiter=",")
        config_list = list(reader)
    return config_list


def run_snippy(snippy_exe, threads, out_dir, ref_genome, r1_seq_file, r2_seq_file):
    """
    This function run snippy software
    :param snippy_exe: the executable of snippy
    :param threads: number of threads
    :param out_dir: the output directory path
    :param ref_genome: the reference genome
    :param r1_seq_file: the R1 fastq file path
    :param r2_seq_file: the R2 fastq file path
    :return: log message
    """
    cmd = "{0} --cpus {1} --outdir {2} --reference {3} --R1 {4} --R2 {5} ".format(snippy_exe, threads, out_dir,
                                                                                  ref_genome,
                                                                                  r1_seq_file, r2_seq_file)
    log_message = " ".join(cmd)
    p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = p.communicate()
    log_message = log_message + '\n' + out.decode("utf-8") + '\n' + err.decode("utf-8")

    print(log_message)

    return log_message


def run_snippy_core_custom(snippy_exe, ref_genome, prefix, bed_file, snippy_folder):
    """
    This function run snippy core custom
    :param snippy_exe: the executable of snippy
    :param ref_genome: the reference genome
    :param prefix: teh prefix given
    :param snippy_folder: snippy folder
    :return: nothing
    """
    cmd = '{0} --ref {1} --prefix {2} --mask {3} {4}'.format(snippy_exe, ref_genome, prefix, bed_file, snippy_folder)
    os.system(cmd)


def get_snippy_dir(geno_ref_dir, result_dir, config_list):
    """
    This function get hash table of snippy folder
    :param geno_ref_dir: the reference genome directory
    :param result_dir: the directory containing the snippy folder
    :param config_list: the configuration list
    :return: hash table of snippy directory
    """
    snippy_dir_dict = {}

    for row in config_list:
        genome_name = row["genomes"].split(".")[0]
        out_dir_root = os.path.join(result_dir, genome_name)
        list_file = os.listdir(out_dir_root)

        out_dir = ""

        for file in list_file:
            file_path = os.path.join(out_dir_root, file)

            if "_" in file:

                file_split = file.split("_")
                for file_part in file_split:
                    if row["strains"] == file_part and os.path.isdir(file_path):
                        out_dir = file_path
                        break
            else:
                if row["strains"] == file and os.path.isdir(file_path):
                    out_dir = file_path
                    break

        if not os.path.exists(os.path.dirname(out_dir)):
            print("Path : {0}".format(os.path.dirname(out_dir)))
            print("ERROR: the directory snippy for the strain {0} dont exist ! Exit!".format(row["strains"]))
            exit(1)
        ref_genome = os.path.join(geno_ref_dir, row["genomes"])

        if ref_genome not in snippy_dir_dict:
            snippy_dir_dict[ref_genome] = [{"out_dir": out_dir, "strain": row["strains"]}]
        else:
            value_list = snippy_dir_dict[ref_genome]

            value_list = value_list + [{"out_dir": out_dir, "strain": row["strains"]}]

            # update dict
            snippy_dir_dict[ref_genome] = value_list

    return snippy_dir_dict


def manage_snippy(read_dir, geno_ref_dir, result_dir, config_list):
    """
    This function manage run of snippy
    :param read_dir: the reads directory
    :param geno_ref_dir: the genome reference directory
    :param result_dir: the directory containing the snippy folder
    :param config_list: the configuration list
    :return: hash table of snippy directory
    """
    print("\n")

    snippy_exe = "snippy"
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
            if row['strains'] + "_" in strain:
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
            snippy_dir_dict[ref_genome] = [{"out_dir": out_dir, "strain": row["strains"]}]
        else:
            value_list = snippy_dir_dict[ref_genome]

            value_list = value_list + [{"out_dir": out_dir, "strain": row["strains"]}]

            # update dict
            snippy_dir_dict[ref_genome] = value_list

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


def manage_snippy_core(snippy_dir_dict, core_genome_path, bed_file):
    """
    This function run snippy-core
    :param snippy_dir_dict: the hash of snippy directory
    :return: a vcf file list
    """

    snippy_exe = "snippy-core"
    name_dir = "core_genome"

    vcf_list = []
    print("\n")

    for genome_ref, snippy_dir_list in snippy_dir_dict.items():

        snippy_dirs = ""
        for element in snippy_dir_list:
            snippy_dirs = snippy_dirs + " " + element["out_dir"]

        # case with only one sample that is not respectable !
        if len(snippy_dir_list) == 1:
            continue

        prefix_snippy_core = os.path.join(core_genome_path, name_dir)
        run_snippy_core_custom(snippy_exe, genome_ref, prefix_snippy_core, bed_file, snippy_dirs)

    for file in os.listdir(core_genome_path):
        if "{0}.vcf".format(name_dir) in file:

            # rename header like the config file want
            vcf_path = os.path.join(core_genome_path, file)

            with open(vcf_path, "r") as vcf:

                content = vcf.readlines()

                for i, line in enumerate(content):
                    if "#CHROM" in line:

                        for genome_ref, snippy_dir_list in snippy_dir_dict.items():
                            for element in snippy_dir_list:
                                line = line.replace(os.path.basename(element["out_dir"]), element["strain"])
                        content[i] = line

            with open(vcf_path, "w") as vcf_write:
                for line in content:
                    vcf_write.write(line)

            vcf_list.append(vcf_path)

    return vcf_list


def manage_filter_snp(vcf_list):
    """
    This function launch the filtering of the vcf files
    :param vcf_list: a list of vcf file
    :return: a list of vcf filtered file
    """
    min_dist = 20
    out_prefix = "filter"
    filter_keep_vcf_list = []

    for vcf_core_file in vcf_list:

        if not os.path.exists(os.path.join(os.path.dirname(vcf_core_file),
                                           out_prefix + "_" +
                                           str(os.path.basename(vcf_core_file).split(".vcf")[0]) +
                                           "_density_filtered_keep.vcf")):

            filter_SNP_density.main(min_dist, vcf_core_file, out_prefix)

        else:
            print("the filtration is already done for {0}".format(vcf_core_file))

        for file in os.listdir(os.path.dirname(vcf_core_file)):

            if "_filtered_keep.vcf" in file:
                filter_keep_vcf_list.append(os.path.join(os.path.dirname(vcf_core_file), file))

    return filter_keep_vcf_list


def read_low_coverage(snippy_dir_dict, snippy_core_genome_folder):
    low_cov_file_list = []

    merge_bed_file = os.path.join(snippy_core_genome_folder, "merged_bed_file.bed")
    merge_bed_sort_file = os.path.join(snippy_core_genome_folder, "merged_bed_sort_file.bed")
    for genome_ref, strain_snippy_list in snippy_dir_dict.items():
        for element in strain_snippy_list:
            for file in os.listdir(element["out_dir"]):
                if "low_coverage_region_" in file:
                    low_coverage_file = os.path.join(element["out_dir"], file)
                    low_cov_file_list.append(low_coverage_file)

    cmd = "cat {0} > {1} | bedtools sort -i {1} | bedtools merge -i stdin > {2}".format(" ".join(low_cov_file_list), merge_bed_file, merge_bed_sort_file)
    log_message = cmd
    p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = p.communicate()
    log_message = log_message + '\n' + out.decode("utf-8") + '\n' + err.decode("utf-8")

    print(log_message)

    awk_cmd = "awk -F'\t' 'BEGIN{SUM=0}{SUM+=$3-$2 }END{print SUM}'"
    cmd = "cat {0} | {1}".format(merge_bed_sort_file, awk_cmd)

    log_message = cmd
    p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = p.communicate()
    log_message = log_message + '\nCOUNT:' + out.decode("utf-8") + '\n' + err.decode("utf-8")

    print(log_message)

    return merge_bed_sort_file


def manage_annotate_filter_snp(filter_keep_vcf_list, snippy_dir_dict):
    """
    This function launch the annotation of filtered vcf file produced before
    :param filter_keep_vcf_list: a list of vcf file
    :param snippy_dir_dict: a dict of snippy strain output
    :return: nothing
    """
    vcf_strain_folder_list = []
    for genome_ref, vcf_snippy_list in snippy_dir_dict.items():
        for strain_dict in vcf_snippy_list:
            vcf_strain_folder_list.append(strain_dict["out_dir"])

    for vcf_core_file in filter_keep_vcf_list:
        output_file = os.path.join(os.path.dirname(vcf_core_file),
                                   "annotated_{0}".format(os.path.basename(vcf_core_file)))
        if not os.path.exists(os.path.join(output_file)):
            annotate_vcf_snippy_core.main(vcf_core_file, vcf_strain_folder_list, output_file)
        else:
            print("the annotation is already done for {0}".format(vcf_core_file))


def manage_r_matrix(filter_keep_vcf_list):
    """
    This function launch the matrix calculation
    :param filter_keep_vcf_list: a list of vcf filtered
    :return: a list of matrix file
    """
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


def manage_make_tree(r_matrix_list, config_list):
    """
    This function create tree
    :param r_matrix_list: a list of matrix file
    :param config_list: the configuration data
    :return: message
    """
    for r_matrix in r_matrix_list:

        filename_newick = "newick_tree.nwk"
        file_newick_path = os.path.join(os.path.dirname(r_matrix), filename_newick)

        if not os.path.exists(file_newick_path):

            with open(r_matrix, "r") as conf:
                reader = DictReader(conf, delimiter="\t")
                headers = reader.fieldnames

                # remove empty col
                headers = list(filter(None, headers))

                matrix = []
                matrix_dict = {}
                for row in reader:

                    line = []
                    for col in headers:
                        value = row[col]
                        line.append(value)
                    matrix.append(line)

                    key = row[""]
                    del row[""]
                    matrix_dict[key] = row

                dm = DistanceMatrix(matrix, headers)

                tree = nj(dm)
                # print(tree.ascii_art())
                # rooted = tree.root_at_midpoint()
                # print(rooted)
                # print(rooted.ascii_art())

                """
                print(headers)
                print(matrix)
                matrix_triangular = []
                for n, line in enumerate(matrix):
                    print(n)
                    desired_array = [int(numeric_string) for numeric_string in line]
                    print(desired_array[0:n+1])
                    matrix_triangular.append(desired_array[0:n+1])

                matrice_biopython = DistanceMatrix(names=headers, matrix=matrix_triangular)
                constructor = DistanceTreeConstructor()
                upgmatree = constructor.upgma(matrice_biopython)

                upgmatree_rooted_at_midpoint = upgmatree.root_at_midpoint()

                newick_str = upgmatree_rooted_at_midpoint
                
                """

                newick_str = tree.root_at_midpoint()

                with open(file_newick_path, 'w') as out:
                    out.write("{0}".format(newick_str))

                print("the newick generation step is done for {0}".format(r_matrix))

                # make phyloxml file
                file_phylo_xml_path = os.path.join(os.path.dirname(r_matrix), "phyloxml.xml")

                Phylo.convert(file_newick_path, 'newick', file_phylo_xml_path, 'phyloxml')
                print("the phyloxml generation step is done for {0}".format(r_matrix))

                # make extended phyloxml file for phyd3
                file_ext_phyloxml_path = os.path.join(os.path.dirname(r_matrix), "extended_phyloxml.xml")

                get_phyloxml_extended(file_ext_phyloxml_path, file_phylo_xml_path, config_list, matrix_dict)
                print("the extended phyloxml generation step is done for {0}".format(r_matrix))
        else:
            print("the newick generation step is already done for {0}".format(r_matrix))


def get_phyloxml_extended(file_ext_phyloxml_path, file_phyloxml_path, config_list, matrix_dict):
    """
    This function make a phylo extended file
    :param file_ext_phyloxml_path: the path of the phylo extended file
    :param file_phyloxml_path: the path of the phylo file
    :param config_list: the configuration data
    :param matrix_dict: the matrix in hash format
    :return: nothing
    """
    strain_label_list = []
    strain_list = []

    pivot_mlst_1 = False
    pivot_mlst_2 = False
    pivot_date_sample = False
    pivot_location = False

    for dict_config in config_list:

        if dict_config.get("MLST-1"):
            pivot_mlst_1 = True
        if dict_config.get("MLST-2"):
            pivot_mlst_2 = True
        if dict_config.get("Date sample"):
            pivot_date_sample = True
        if dict_config.get("Location"):
            pivot_location = True

    with open(file_phyloxml_path, 'r') as xml_very_uggly:
        content = xml_very_uggly.readlines()
        for n, line in enumerate(content):
            newline = line.lstrip().rstrip()
            content[n] = newline
            if "<name>" in newline and "root" not in newline:
                strain = newline.replace("<name>", "").replace("</name>", "")
                strain_list.append(strain)
                if strain[0].isnumeric():
                    strain = "s" + strain
                strain_label_list.append(strain)

    xml_string = "".join(content)

    tree = ET.ElementTree(ET.fromstring(xml_string))
    root = tree.getroot()

    ################################################
    # adding Label

    # adding an element to the root node
    labels = ET.SubElement(root, "labels")

    # adding an element to the labels node
    attrib = {'type': 'text'}

    # ST 1
    if pivot_mlst_1:
        label = ET.SubElement(labels, "label", attrib)
        ET.SubElement(label, "name", show="1").text = "ST-1"
        ET.SubElement(label, "data", tag="ST-1")

    # ST 2
    if pivot_mlst_2:
        label = ET.SubElement(labels, "label", attrib)
        ET.SubElement(label, "name", show="1").text = "ST-2"
        ET.SubElement(label, "data", tag="ST-2")

    # Date sample
    if pivot_date_sample:
        label = ET.SubElement(labels, "label", attrib)
        ET.SubElement(label, "name", show="1").text = "Date sample"
        ET.SubElement(label, "data", tag="Date-sample")

    # Location
    if pivot_location:
        label = ET.SubElement(labels, "label", attrib)
        ET.SubElement(label, "name", show="1").text = "location"
        ET.SubElement(label, "data", tag="location")

    for strain in strain_label_list:
        label = ET.SubElement(labels, "label", attrib)
        ET.SubElement(label, "name", show="1").text = strain
        ET.SubElement(label, "data", tag=strain)

    #######################################
    # add labels

    leaf = './{http://www.phyloxml.org}phylogeny/{http://www.phyloxml.org}clade'

    resu_dict = get_dict_strain_et(root, leaf, {}, [])

    for leaf, strain in resu_dict.items():
        for config_dict in config_list:
            if config_dict.get("strains") == strain:

                if pivot_mlst_1:
                    ET.SubElement(leaf, "{http://www.phyloxml.org}ST-1").text = config_dict.get("MLST-1")
                if pivot_mlst_2:
                    ET.SubElement(leaf, "{http://www.phyloxml.org}ST-2").text = config_dict.get("MLST-2")
                if pivot_date_sample:
                    ET.SubElement(leaf, "{http://www.phyloxml.org}Date-sample").text = config_dict.get("Date sample")
                if pivot_location:
                    ET.SubElement(leaf, "{http://www.phyloxml.org}location").text = config_dict.get("Location")

                for n, strain_label in enumerate(strain_label_list):
                    ET.SubElement(leaf, ("{http://www.phyloxml.org}"+strain_label)).text =\
                        str(matrix_dict.get(strain).get(strain_list[n]))

                break
            else:
                continue

    ################################################
    tree.write(file_ext_phyloxml_path, encoding="utf-8")
    ############
    # make xml pretty

    xml_str = minidom.parse(file_ext_phyloxml_path).toprettyxml(indent="   ")
    xml_str = xml_str.replace('phy:', '')
    
    with open(file_ext_phyloxml_path, "w") as f:
        f.write(xml_str)


def get_dict_strain_et(root, leaf, resu_dict, leaf_done_list):
    """
    This function render a dict of a ET tree
    :param root: teh root of the ET tree
    :param leaf: the first leaf
    :param resu_dict: the result dict
    :param leaf_done_list: the list of leaf done in tree
    :return: the result dict
    """

    if leaf not in leaf_done_list:
        if root.findall(leaf):
            for clade in root.findall(leaf):
                name = clade.find('{http://www.phyloxml.org}name')

                if name is not None and clade.find('{http://www.phyloxml.org}name').text == "root":
                    resu_dict = get_dict_strain_et(root, leaf + "/" + clade.tag, resu_dict, leaf_done_list)
                elif name is not None:
                    resu_dict[clade] = clade.find('{http://www.phyloxml.org}name').text
                    leaf_done_list.append(leaf)
                else:
                    resu_dict = get_dict_strain_et(root, leaf + "/" + clade.tag, resu_dict, leaf_done_list)
    return resu_dict


def manage_snp_network(r_matrix_list, config_file, filter_keep_vcf_list):
    """
    This function make the networkx/bokeh tree
    :param r_matrix_list: the list of matrix
    :param config_file: the configuration file
    :param filter_keep_vcf_list: the list of filtered vcf files
    :return: nothing
    """
    for mtx_file in r_matrix_list:
        graph_name = os.path.join(os.path.dirname(mtx_file), "SNP_network.html")

        with open(filter_keep_vcf_list[0]) as f:
            count_snp_keep = sum(1 for line in f) - 6

        mtx2mst.main(mtx_file, graph_name, config_file, count_snp_keep)


def pre_main(args):
    """
    The pre-main function that recives arguments
    :param args: the arguments
    :return: the arguments
    """
    read_dir = geno_ref_dir = result_dir = config_file = ""
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
    :param jump_snippy_detection:
    :param read_dir:
    :param geno_ref_dir:
    :param result_dir:
    :param config_file:
    """
    config_list = load_config(config_file)

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

    name_dir = "core_genome"
    for geno_ref, snippy_dir_list in snippy_dir_dict.items():
        snippy_core_genome_folder = os.path.join(os.path.dirname(snippy_dir_list[0]["out_dir"]), name_dir)
        break

    if not os.path.exists(snippy_core_genome_folder):
        os.makedirs(snippy_core_genome_folder)

    # filter by low coverage
    print("\nStart load filter by low coverage")
    bed_file = read_low_coverage(snippy_dir_dict, snippy_core_genome_folder)
    print("\nEnd load filter by low coverage")
    print("*********************************************")

    # snippy_core
    print("\nStart Snippy-core")
    vcf_list = manage_snippy_core(snippy_dir_dict, snippy_core_genome_folder, bed_file)
    print("End Snippy-core")
    print("*********************************************")

    # filter SNP
    print("\nStart filtering SNP")
    filter_keep_vcf_list = manage_filter_snp(vcf_list)
    print("End filtering SNP")
    print("*********************************************")

    # annotated filter SNP
    print("\nStart annotate filtering SNP")
    manage_annotate_filter_snp(filter_keep_vcf_list, snippy_dir_dict)
    print("End annotate filtering SNP")
    print("*********************************************")

    # R_matrix
    print("\nStart distance matrix")
    r_matrix_list = manage_r_matrix(filter_keep_vcf_list)
    print("End distance matrix")
    print("*********************************************")

    # Make tree
    print("\nStart make tree")
    manage_make_tree(r_matrix_list, config_list)
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

    usage = "CNR Phylo tree is a snakemake based program"

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
    """
    The version of the script
    :return: the version of the script
    """
    return "1.0"


if __name__ == '__main__':
    run()
