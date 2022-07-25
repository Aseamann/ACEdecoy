import argparse
import subprocess
import time
import os
from shutil import copyfile


rosetta_dir = ""
program_dir = os.path.dirname(os.path.realpath(__file__))
version = ""
relax_protein_runs = 100
relax_xml_runs = 40
relax_bb_runs = 30
relax_fast_runs = 30
docking_runs = 5000
refine_runs = 100

cpu_protein_relax = 51
cpu_xml_relax = 41
cpu_bb_relax = 31
cpu_fast_relax = 31
cpu_dock = 81


# Controls the executions of rosetta
def run_rosetta(pdb):
    print("Starting!")
    print("Running relax...")
    start = time.time()
    make_dirs()
    run_protein_relax(pdb)
    print("Best Refinement: " + check_score_protein_relax())
    print("DONE!")
    print(f"Time: {(time.time() - start)/60:.0f} mins")


def make_dirs():
    os.mkdir(program_dir + "/output_files/")
    os.mkdir(program_dir + "/output_files/relax/")


###############################################
# pMHC
# Runs protein relax
def run_protein_relax(protein):
    dir_relax = rosetta_dir + "/main/source/bin/relax.mpi." + version
    make_protein_relax_file(protein)
    subprocess.run(["mpirun", dir_relax, "@flag_protein_relax"], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    check_score_protein_relax()


# Checks for best protein relaxed structure. Returns name and copies it from output to input.
def check_score_protein_relax():
    score_dic = {}
    best_pdb = ""
    score = 100000.0
    with open("output_files/relax/score.sc", "r") as score_read:
        for line in score_read:
            if not line.__contains__("SEQUENCE") and not line.__contains__("total_score"):
                score_dic[line.split()[-1]] = line.split()
    for pdb in score_dic:
        if float(score_dic[pdb][1]) < score:
            best_pdb = score_dic[pdb][-1]
            score = float(score_dic[pdb][1])
    copyfile("output_files/relax/" + best_pdb + ".pdb", best_pdb + ".pdb")
    return best_pdb


def make_protein_relax_file(protein):
    with open("flag_protein_relax", "w") as relax_file:
        relax_file.write("-in:file:s " + protein + "\n\n")
        relax_file.write("#SBATCH --ntasks=" + str(cpu_protein_relax) + "\n")
        relax_file.write("-nstruct " + str(relax_protein_runs) + " \n\n")
        relax_file.write("#-relax:constrain_relax_to_start_coords\n")
        relax_file.write("#-relax:ramp_constraints false\n\n")
        relax_file.write("-ex1\n-ex2\n\n")
        relax_file.write("-use_input_sc\n-flip_HNQ\n-no_optH false\n\n")
        relax_file.write("-out:path:all output_files/relax\n")


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("pdb", help="PDB file", type=str)
    parser.add_argument("-l", "--linux", help="Changes to linux runnable program.", default=True,
                        dest='linux', action='store_true')
    parser.add_argument("-m", "--mac", help="Changes to mac runnable program.", default=False,
                        dest='mac', action='store_true')
    return parser.parse_args()


def main():
    args = parse_args()
    global version
    if args.mac:
        version = "macosclangrelease"
    elif args.linux:
        version = "linuxgccrelease"
    global rosetta_dir
    rosetta_dir = "/Users/austinseamann/Rosetta/rosetta.source.release-314/"
    run_rosetta(args.pdb)


if __name__ == '__main__':
    main()

