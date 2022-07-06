import shutil
import subprocess
from datetime import datetime


def main(file_newick_path, file_png_path):

    # R Script
    print(f'Execute RScript {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}')
    ex_r = shutil.which("Rscript")
    r_script = "/usr/local/CNRPhyloTree/cnr_phylo_tree_src/save_newick_to_png.R"
    print(r_script)

    cmd = f"{ex_r} {r_script} -i {file_newick_path} -o {file_png_path}"
    log_message = f"Command used : \n {cmd}\n"
    # launch
    process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT).stdout.read()
    log_message = log_message + str(process)
    print(log_message)
