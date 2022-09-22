#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
    This script will manage the phylo tree construction
"""

import argparse
import csv
import multiprocessing
import os
import subprocess
import sys
from csv import DictReader
from xml.dom import minidom

from Bio import Phylo, SeqIO
# from Bio.Phylo.TreeConstruction import DistanceTreeConstructor, DistanceMatrix
from skbio import DistanceMatrix
from skbio.tree import nj

from cnr_phylo_tree_src import filter_SNP_density, vcf2dist, mtx2mst, annotate_vcf_snippy_core, pre_filter_SNP_density, \
    execute_R
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
    cmd = f"{snippy_exe} --cpus {threads} --outdir {out_dir} --reference {ref_genome} --R1 {r1_seq_file}" \
          f" --R2 {r2_seq_file} "
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
    if bed_file == "no file":
        cmd = f'{snippy_exe} --ref {ref_genome} --prefix {prefix} {snippy_folder}'
    else:
        cmd = f'{snippy_exe} --ref {ref_genome} --prefix {prefix} --mask {bed_file} {snippy_folder}'
    print(cmd)
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
    ref_genome = ""
    for row in config_list:
        genome_name = row["genomes"].split(".")[0]
        out_dir_root = os.path.join(result_dir, genome_name)
        list_file = os.listdir(out_dir_root)
        out_dir = ""
        for file in list_file:
            file_path = os.path.join(out_dir_root, file)
            if "_" in str(file):
                file_split = str(file).split("_")
                for file_part in file_split:
                    if row["strains"] == file_part and os.path.isdir(file_path):
                        out_dir = file_path
                        break
            else:
                if row["strains"] == file and os.path.isdir(file_path):
                    out_dir = file_path
                    break
        if not os.path.exists(os.path.dirname(out_dir)):
            print(f"Path : {out_dir}")
            print(f"ERROR: the directory snippy for the strain {row['strains']} dont exist ! Exit!")
            exit(1)
        # case if reference genome in analysis
        ref_genome = os.path.join(geno_ref_dir, row["genomes"])
        print(ref_genome)
        if ref_genome not in snippy_dir_dict:
            snippy_dir_dict[ref_genome] = [{"out_dir": out_dir, "strain": row["strains"]}]
        else:
            """
            value_list = snippy_dir_dict[ref_genome]
            value_list = value_list + [{"out_dir": out_dir, "strain": row["strains"]}]
            # update dict
            snippy_dir_dict[ref_genome] = value_list
            """
            print("lol")
    print(snippy_dir_dict, ref_genome)
    return snippy_dir_dict, ref_genome


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
        r1_seq_file = r2_seq_file = ""
        ref_genome = os.path.join(geno_ref_dir, row["genomes"])
        if not os.path.exists(ref_genome):
            print(f"your reference genome file {ref_genome} don't exists or we haven't the permission\n")
            print(usage)
            sys.exit()
        for strain in strain_list:
            if row['strains'] + "_" in strain:
                if "R1" in strain:
                    r1_seq_file = os.path.join(read_dir, strain)
                    if not os.path.exists(r1_seq_file):
                        print(f"the fastq file {r1_seq_file} don't exists or we haven't the permission\n")
                        print(usage)
                        sys.exit()
                if "R2" in strain:
                    r2_seq_file = os.path.join(read_dir, strain)
                    if not os.path.exists(r2_seq_file):
                        print(f"the fastq file {r2_seq_file} don't exists or we haven't the permission\n")
                        print(usage)
                        sys.exit()
        if not r1_seq_file:
            print(f"The sequence file {r1_seq_file} is not found")
            print(f"The program continue but without the strain  {row['strains']}")
            continue
        if not r2_seq_file:
            print(f"The sequence file {r1_seq_file} is not found")
            print(f"The program continue but without the strain  {row['strains']}")
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
                name=f'snippy {row["strains"]}'
            )
            jobs.append(p)
        else:
            print(f"The snippy folder {out_dir} already exist !")
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
    return snippy_dir_dict, ref_genome


def manage_snippy_core(snippy_dir_dict, core_genome_path, bed_file):
    """
    This function run snippy-core
    :param snippy_dir_dict: the hash of snippy directory
    :return: a vcf file
    """
    snippy_exe = "snippy-core"
    name_dir = "core_genome"
    vcf_path = ""
    print("\n")
    print("snippy_dir_dict : ", snippy_dir_dict)
    print("core_genome_path : ", core_genome_path)
    print("Bed file : ", bed_file)
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
        # rename header like the config file want
        if f"{name_dir}.vcf" not in file:
            continue
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
    return vcf_path


def manage_filter_snp(vcf_path, min_dist, type_matrice):
    """
    This function launch the filtering of the vcf files
    :return: a list of vcf filtered file
    """
    out_prefix = "filter"
    filter_keep_vcf_list = []
    bool = False
    if not os.path.exists(os.path.join(os.path.dirname(vcf_path), out_prefix + "_" + str(os.path.basename(vcf_path).split(".vcf")[0]) + "_density_filtered_keep.vcf")):
        bool = filter_SNP_density.main(min_dist, vcf_path, out_prefix, type_matrice)
    else:
        print(f"the filtration is already done for {vcf_path}")
    if bool:
        for file in os.listdir(os.path.dirname(vcf_path)):
            if "_filtered_keep.vcf" in file:
                filter_keep_vcf_list.append(os.path.join(os.path.dirname(vcf_path), file))
    return filter_keep_vcf_list


def manage_pre_filter_snp(vcf_path):
    """
    This function launch the pre-filtering of the vcf files
    :return: a list of vcf filtered file
    """
    out_prefix = "pre-filter"
    dir_name_path = os.path.dirname(vcf_path)
    name_file_vcf_input = os.path.basename(os.path.splitext(vcf_path)[0])
    name_output_file = out_prefix + "_" + name_file_vcf_input
    out_vcf_file = os.path.join(dir_name_path, name_output_file + "_pre-keep.vcf")
    out_vcf_unkeep_file = os.path.join(dir_name_path, name_output_file + "_pre-unkeep.vcf")
    pre_filter_SNP_density.main(vcf_path, out_vcf_file, out_vcf_unkeep_file)
    return out_vcf_file


def read_low_coverage(snippy_dir_dict, snippy_core_genome_folder):
    low_cov_file_list = []
    merge_bed_file = os.path.join(snippy_core_genome_folder, "merged_bed_file.bed")
    merge_bed_sort_file = os.path.join(snippy_core_genome_folder, "merged_bed_sort_file.bed")
    for genome_ref, strain_snippy_list in snippy_dir_dict.items():
        for element in strain_snippy_list:
            low_coverage_file = ""
            low_coverage_file_merge = ""
            for file in os.listdir(element["out_dir"]):
                if "low_coverage_region_" in file:
                    low_coverage_file = os.path.join(element["out_dir"], file)
                elif "merge_bad_position_sort" in file:
                    low_coverage_file_merge = os.path.join(element["out_dir"], file)
            if low_coverage_file_merge:
                low_cov_file_list.append(low_coverage_file_merge)
            elif low_coverage_file:
                low_cov_file_list.append(low_coverage_file)

    if low_cov_file_list:
        print("LOW COVERAGE PROCESS")
        # cat
        print("Length files list", len(low_cov_file_list))
        x = 0
        tmp_cat_list = []
        for i in range(0, len(low_cov_file_list), 50):
            tmp_cat_file = os.path.join(snippy_core_genome_folder, f"merged_bed_file_{x}.bed")
            cmd = f"cat {' '.join(low_cov_file_list[i:i + 50])} > {tmp_cat_file} "
            log_message = cmd
            p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            out, err = p.communicate()
            log_message = log_message + '\n' + out.decode("utf-8") + '\n' + err.decode("utf-8")
            print(log_message)
            tmp_cat_list.append(tmp_cat_file)
            x += 1

        cmd = f"cat {' '.join(tmp_cat_list)} > {merge_bed_file} "
        log_message = cmd
        p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        out, err = p.communicate()
        log_message = log_message + '\n' + out.decode("utf-8") + '\n' + err.decode("utf-8")
        print(log_message)
        # bedtools
        cmd = f"bedtools sort -i {merge_bed_file} | bedtools merge -i stdin > {merge_bed_sort_file}"
        log_message = cmd
        p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        out, err = p.communicate()
        log_message = log_message + '\n' + out.decode("utf-8") + '\n' + err.decode("utf-8")
        print(log_message)
        awk_cmd = "awk -F'\t' 'BEGIN{SUM=0}{SUM+=$3-$2 }END{print SUM}'"
        cmd = f"cat {merge_bed_sort_file} | {awk_cmd}"
        log_message = cmd
        p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        out, err = p.communicate()
        log_message = log_message + '\nCOUNT:' + out.decode("utf-8") + '\n' + err.decode("utf-8")
        print(log_message)
        return merge_bed_sort_file
    else:
        return "no file"


def compute_constant_site(bed_file, ref_genome):
    ref_A = ref_T = ref_G= ref_C= ref_O = constant_A = constant_T =constant_G =constant_C =constant_O = 0
    if bed_file != "no file":
        out_constant_fasta = os.path.join(os.path.dirname(bed_file), "constant_sites.fasta")
        if ".gbff" in os.path.basename(ref_genome) or ".gbk" in os.path.basename(ref_genome):
            new_ref_genome = os.path.join(os.path.dirname(ref_genome), os.path.basename(ref_genome).split(".")[0]+".fasta")
            cmd = f"any2fasta {ref_genome} > {new_ref_genome}"
            log_message = cmd
            p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            out, err = p.communicate()
            log_message = log_message + '\n' + out.decode("utf-8") + '\n' + err.decode("utf-8")
            print(log_message)
            ref_genome = new_ref_genome

        for rec in SeqIO.parse(ref_genome, "fasta"):
            for base in rec.seq:
                if base.upper() == "A":
                    ref_A += 1
                elif base.upper() == "T":
                    ref_T += 1
                elif base.upper() == "G":
                    ref_G += 1
                elif base.upper() == "C":
                    ref_C += 1
                else:
                    ref_O += 1

        cmd = f"bedtools getfasta -fi {ref_genome} -bed {bed_file} > {out_constant_fasta}"
        log_message = cmd
        p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        out, err = p.communicate()
        log_message = log_message + '\n' + out.decode("utf-8") + '\n' + err.decode("utf-8")
        print(log_message)

        for rec in SeqIO.parse(out_constant_fasta, "fasta"):
            for base in rec.seq:
                if base.upper() == "A":
                    constant_A += 1
                elif base.upper() == "T":
                    constant_T += 1
                elif base.upper() == "G":
                    constant_G += 1
                elif base.upper() == "C":
                    constant_C += 1
                else:
                    constant_O += 1

        out_constant_tsv = os.path.join(os.path.dirname(bed_file), 'constant_site.tsv')
        with open(out_constant_tsv, "w") as out_file:
            writer = csv.writer(out_file, delimiter='\t')
            writer.writerow(["genomes", "A", "T", "G", "C", "sum"])
            sum_ref = ref_A+ref_T+ref_G+ref_C
            sum_constant = sum_ref - (constant_A+constant_T+constant_G+constant_C)
            writer.writerow([os.path.basename(ref_genome).split(".")[0], ref_A, ref_T, ref_G, ref_C, sum_ref])
            writer.writerow([os.path.basename(out_constant_fasta).split(".")[0], ref_A-constant_A, ref_T-constant_T,
                             ref_G-constant_G, ref_C-constant_C, sum_constant])


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
    print("vcf_strain_folder_list: ", vcf_strain_folder_list)
    for vcf_core_file in filter_keep_vcf_list:
        print("vcf_core_file: ", vcf_core_file)
        output_file = os.path.join(os.path.dirname(vcf_core_file), f"annotated_{os.path.basename(vcf_core_file)}")
        if not os.path.exists(os.path.join(output_file)):
            annotate_vcf_snippy_core.main(vcf_core_file, vcf_strain_folder_list, output_file)
        else:
            print(f"the annotation is already done for {vcf_core_file}")


def manage_r_matrix(filter_keep_vcf_list, type_matrice):
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
                print(f"the matrix step is already done for {filter_keep_vcf}")
    if not jump:
        for filter_keep_vcf in filter_keep_vcf_list:
            vcf2dist.main(filter_keep_vcf, type_matrice)
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
        filename_png = "tree.png"
        file_newick_path = os.path.join(os.path.dirname(r_matrix), filename_newick)
        file_png_path = os.path.join(os.path.dirname(r_matrix), filename_png)
        if not os.path.exists(file_newick_path):
            with open(r_matrix, "r") as conf:
                reader = DictReader(conf, delimiter="\t")
                headers = reader.fieldnames
                # remove empty col
                headers = list(filter(None, headers))
                matrix = []
                matrix_dict = {}
                x = 0

                for row in reader:
                    x += 1
                    y = 0
                    line = []
                    for col in headers:
                        y += 1
                        value = row[col]
                        line.append(value)
                    print(x, y)
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
                    out.write(f"{newick_str}")
                print(f"the newick generation step is done for {r_matrix}")
                # make phyloxml file
                file_phylo_xml_path = os.path.join(os.path.dirname(r_matrix), "phyloxml.xml")
                Phylo.convert(file_newick_path, 'newick', file_phylo_xml_path, 'phyloxml')
                print(f"the phyloxml generation step is done for {r_matrix}")
                # make extended phyloxml file for phyd3
                file_ext_phyloxml_path = os.path.join(os.path.dirname(r_matrix), "extended_phyloxml.xml")
                get_phyloxml_extended(file_ext_phyloxml_path, file_phylo_xml_path, config_list, matrix_dict)
                print(f"the extended phyloxml generation step is done for {r_matrix}")

                print()

                execute_R.main(file_newick_path, file_png_path)

        else:
            print(f"the newick generation step is already done for {r_matrix}")


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
    pivot_mlst_1 = pivot_mlst_2 = pivot_date_sample = pivot_location = False
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
        ET.SubElement(label, "name", show="1").text = "Location"
        ET.SubElement(label, "data", tag="Location")
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
                    ET.SubElement(leaf, "{http://www.phyloxml.org}Location").text = config_dict.get("Location")
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
    read_dir = geno_ref_dir = result_dir = config_file = type_matrice = ""
    min_dist = 0
    jump_snippy_detection = filter_homogeneous_SNP = False
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
        print(f"Result directory was created at '{result_dir}'")
        os.makedirs(result_dir)
    if args.FilterHomogeneousSNP:
        filter_homogeneous_SNP = args.FilterHomogeneousSNP
    if args.minimum_distance:
        min_dist = args.minimum_distance
    if args.type_matrice:
        type_matrice = args.type_matrice

    main(read_dir, geno_ref_dir, result_dir, config_file, min_dist, type_matrice, jump_snippy_detection,
         filter_homogeneous_SNP)


def main(read_dir, geno_ref_dir, result_dir, config_file, min_dist, type_matrice, jump_snippy_detection=False,
         filter_homogeneous_SNP=False):
    """
    This make a configuration file with argument for launch of snakemake
    :param jump_snippy_detection:
    :param read_dir:
    :param geno_ref_dir:
    :param result_dir:
    :param config_file:
    """

    print("#################\n")
    print(f"Jump Snippy Detection : {jump_snippy_detection} \n")
    print("#################\n")
    print(f"Filter Homogeneous SNP : {filter_homogeneous_SNP} \n")
    print("#################\n")
    ref_genome = ""
    config_list = load_config(config_file)
    if not jump_snippy_detection:
        # snippy
        print("\nStart Snippy")
        snippy_dir_dict, ref_genome = manage_snippy(read_dir, geno_ref_dir, result_dir, config_list)
        print("End Snippy")
        print("*********************************************")
    else:
        print("Skip Snippy Detection")
        snippy_dir_dict, ref_genome = get_snippy_dir(geno_ref_dir, result_dir, config_list)
        print("END Skip Snippy Detection")
        print("*********************************************")
    name_dir = "core_genome"
    snippy_core_genome_folder = ""
    number_strain = 0
    for geno_ref, snippy_dir_list in snippy_dir_dict.items():
        number_strain = len(snippy_dir_list)
        snippy_core_genome_folder = os.path.join(os.path.dirname(snippy_dir_list[0]["out_dir"]), name_dir)
        break
    if not os.path.exists(snippy_core_genome_folder):
        os.makedirs(snippy_core_genome_folder)
    print(f"Number Strain = {number_strain}")
    # filter by low coverage
    print("\nStart load filter by low coverage")
    bed_file = read_low_coverage(snippy_dir_dict, snippy_core_genome_folder)
    print("\nEnd load filter by low coverage")
    print("*********************************************")
    # calculate constant site
    print("\nStart compute constant site")
    compute_constant_site(bed_file, ref_genome)
    print("\nEnd compute constant site")
    print("*********************************************")
    # snippy_core
    print("\nStart Snippy-core")
    vcf_path = manage_snippy_core(snippy_dir_dict, snippy_core_genome_folder, bed_file)
    print("End Snippy-core")
    print("*********************************************")
    # pre-filter SNP
    if filter_homogeneous_SNP:
        print("\nStart pre-filtering SNP")
        vcf_path = manage_pre_filter_snp(vcf_path)
        print("End pre-filtering SNP")
        print("*********************************************")
    # filter SNP
    print("\nStart filtering SNP")
    filter_keep_vcf_list = manage_filter_snp(vcf_path, min_dist, type_matrice)
    print("End filtering SNP")
    print("*********************************************")

    if filter_keep_vcf_list:
        # annotated filter SNP
        print("\nStart annotate filtering SNP")
        manage_annotate_filter_snp(filter_keep_vcf_list, snippy_dir_dict)
        print("End annotate filtering SNP")
        print("*********************************************")
        # R_matrix
        print("\nStart distance matrix")
        r_matrix_list = manage_r_matrix(filter_keep_vcf_list, type_matrice)
        print("End distance matrix")
        print("*********************************************")
        if number_strain >= 3:
            # Make tree
            print("\nStart make tree")
            manage_make_tree(r_matrix_list, config_list)
            print("End make tree")
            print("*********************************************")
            # networkX
            print("\nStart networkX")
            manage_snp_network(r_matrix_list, config_file, filter_keep_vcf_list)
            print("End networkX")
        else:
            print('The number of strain to make a tree is too low')
    else:
        print("No SNP are in the -keep- list. The program is stopped.")


def run():
    """
    Take the argument of the command and give a variable for this
    :return: args
    """
    global usage
    usage = "CNR Phylo tree is a snakemake based program"
    parser = argparse.ArgumentParser(
        description="This program make a phylogenetic tree with NGS sequence reads and one or more reference genome")
    parser.add_argument('-i', '--repRead', dest="repRead", default='', help='Enter the path of the directory which'
                                                                            ' contain your fasta/fastq/fas/fa(.gz)'
                                                                            ' reads')
    parser.add_argument('-o', '--repResult', dest="repResult", default='',
                        help="Enter the path of your directory where you want your result (without other result or"
                             " risk of nothing that will be done)")
    parser.add_argument('-g', '--repGenomeRef', dest="genomeRef", default='',
                        help="Enter the path of directory with reference genome (if no configfile enter the path to"
                             " your reference genome)")
    parser.add_argument('-c', '--config', dest="config", default='',
                        help="Enter the path of your csv file which config the association of genome with"
                             "strain (help of csv2json.py for see how make this csv file)")
    parser.add_argument('-md', '--minimum_distance', dest="minimum_distance", default=20,
                        help="The minimum distance between 2 SNP to consider as recombinant.")
    parser.add_argument('-tm', '--type_matrice', dest="type_matrice", default="absolu",
                        help="The mode of computation of the distance matrice. [absolu, relatif, relatif-norm]"
                             " (default: absolu)")
    parser.add_argument('-p', '--SnippyAlreadyDone', dest="Already", action='store_true', default=False,
                        help="Indicate if all snippy detection are already done (default: False)")
    parser.add_argument('-f1', '--FilterHomogeneousSNP', dest="FilterHomogeneousSNP",
                        action='store_true', default=False,
                        help="Allow to filter homogeneous SNP before filtering (default: False)")

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
