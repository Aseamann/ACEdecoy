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
    # run_xml_relax("peptide.pdb")
    # run_bb_relax("peptide.pdb")
    # run_fast_relax("peptide.pdb")
    # print("Preparing prepack...")
    # run_prepack(pdb)
    # print("Running docking...")
    # run_dock(pdb)
    # remove_dock(check_score_dock())
    # print("Running refine...")
    # pdb_refine = check_score_dock()
    # run_refine(pdb_refine + ".pdb")
    # best_refine = check_score_refine()
    # remove_refine(best_refine)
    # print(best_refine)
    print("DONE!")
    print(f"Time: {(time.time() - start)/60:.0f} mins")


def make_dirs():
    # os.mkdir(program_dir + "/input_files/")
    os.mkdir(program_dir + "/output_files/")
    # os.mkdir(program_dir + "/output_files/dock/")
    os.mkdir(program_dir + "/output_files/relax/")
    # os.mkdir(program_dir + "/output_files/refine/")
    # os.mkdir(program_dir + "/output_files/prepack/")
    # os.mkdir(program_dir + "/input_files/protein_ensembles/")
    # os.mkdir(program_dir + "/input_files/peptide_ensembles/")


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


###############################################
# XML
# Runs xml relax for peptide
def run_xml_relax(peptide):
    dir_script = rosetta_dir + "/main/source/bin/rosetta_scripts.mpi." + version
    make_xml_file()
    make_xml_flag(peptide)
    subprocess.run(["mpirun", dir_script, "@flag_xml_relax"], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)


# Creates the xml file for normal relax run of peptide
def make_xml_file():
    with open("nma.xml", "w") as f:
        f.write("<ROSETTASCRIPTS>\n")
        f.write("\t<SCOREFXNS>\n\t\t<ScoreFunction name=\"bn15_cart\" weights=\"beta_nov15_cart\" />\n")
        f.write("\t</SCOREFXNS>\n\t<RESIDUE_SELECTORS>\n\t</RESIDUE_SELECTORS>\n")
        f.write("\t<TASKOPERATIONS>\n\t</TASKOPERATIONS>\n\t<FILTERS>\n\t</FILTERS>\n")
        f.write("\t<MOVERS>\n\t\t<NormalModeRelax name=\"nma\" cartesian=\"true\" centroid=\"false\"\n")
        f.write("\t\t\tscorefxn=\"bn15_cart\" nmodes=\"5\" mix_modes=\"true\" pertscale=\"1.0\"\n")
        f.write("\t\t\trandomselect=\"false\" relaxmode=\"relax\" nsample=\"120\"\n")
        f.write("\t\t\tcartesian_minimize=\"false\" />\n\t</MOVERS>\n")
        f.write("\t<APPLY_TO_POSE>\n\t</APPLY_TO_POSE>\n")
        f.write("\t<PROTOCOLS>\n\t\t<Add mover=\"nma\" />\n\t</PROTOCOLS>\n")
        f.write("\t<OUTPUT scorefxn=\"bn15_cart\" />\n")
        f.write("</ROSETTASCRIPTS>\n")


# Creates the flag file for the normal xml relax
def make_xml_flag(peptide):
    with open("flag_xml_relax", "w") as f:
        f.write("-in:file:s " + str(peptide) + "\n\n")
        f.write("#SBATCH --ntasks=" + str(cpu_xml_relax) + "\n")
        f.write("-nstruct " + str(relax_xml_runs) + "\n\n")
        f.write("-parser:protocol nma.xml\n\n")
        f.write("-out:path:all input_files/peptide_ensembles\n")
        f.write("-out:suffix _xml\n")


###############################################
# BB Relax
# Runs bb relax for peptide
def run_bb_relax(peptide):
    dir_relax = rosetta_dir + "/main/source/bin/relax.mpi." + version
    make_bb_relax_file(peptide)
    subprocess.run(["mpirun", dir_relax, "@flag_bb_peptide_relax"], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)


def make_bb_relax_file(peptide):
    with open("flag_bb_peptide_relax", "w") as relax_file:
        relax_file.write("-in:file:s " + str(peptide) + "\n\n")
        relax_file.write("#SBATCH --ntasks=" + str(cpu_bb_relax) + "\n")
        relax_file.write("-nstruct " + str(relax_bb_runs) + " \n\n")
        relax_file.write("-backrub:ntrials 20000\n")
        relax_file.write("-backrub:mc_kt 0.6\n\n")
        relax_file.write("-out:path:all input_files/peptide_ensembles/\n")
        relax_file.write("-out:suffix _bb\n")


###############################################
# Fast Relax
# Runs fast relax for peptide
def run_fast_relax(peptide):
    dir_relax = rosetta_dir + "/main/source/bin/relax.mpi." + version
    make_fast_relax_file(peptide)
    subprocess.run(["mpirun", dir_relax, "@flag_fast_relax"], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)


def make_fast_relax_file(peptide):
    with open("flag_fast_relax", "w") as relax_file:
        relax_file.write("-in:file:s " + str(peptide) + "\n\n")
        relax_file.write("#SBATCH --ntasks=" + str(cpu_fast_relax) + "\n")
        relax_file.write("-nstruct " + str(relax_fast_runs) + " \n\n")
        relax_file.write("-relax:thorough\n\n")
        relax_file.write("-out:path:all input_files/peptide_ensembles/\n")
        relax_file.write("-out:suffix _fast\n")


###############################################
# PREPACK
def run_prepack(pdb):
    dir_relax = rosetta_dir + "/main/source/bin/docking_prepack_protocol.mpi." + version
    make_ensemble_files()
    make_prepack_file(pdb)
    subprocess.run(["mpirun", dir_relax, "@flag_ensemble_prepack"], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)


def make_ensemble_files():
    global program_dir
    with open("protein_ensemblelist", "w") as p:
        for filename in sorted(os.listdir(program_dir + "/input_files/protein_ensembles/")):
            if filename.endswith(".pdb"):
                p.write("input_files/protein_ensembles/" + filename + "\n")
    with open("peptide_ensemblelist", "w") as t:
        for filename in sorted(os.listdir(program_dir + "/input_files/peptide_ensembles")):
            if filename.endswith(".pdb"):
                t.write("input_files/peptide_ensembles/" + filename + "\n")


def make_prepack_file(pdb):
    with open("flag_ensemble_prepack", "w") as f:
        f.write("-in:file:s " + str(pdb) + "\n")
        f.write("-unboundrot " + str(pdb) + "\n\n")
        f.write("-nstruct 1\n")
        f.write("-partners E_A\n\n")
        f.write("-ensemble1 protein_ensemblelist\n")
        f.write("-ensemble2 peptide_ensemblelist\n")
        f.write("-ex1\n-ex2aro\n\n")
        f.write("-out:path:all output_files/prepack\n")
        f.write("-out:suffix _prepack\n")


###############################################
# DOCK
def run_dock(pdb):
    dir_dock = rosetta_dir + "/main/source/bin/docking_protocol.mpi." + version
    make_docking_file(pdb)
    subprocess.run(["mpirun", dir_dock, "@flag_ensemble_docking"], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)


def make_docking_file(pdb):
    with open("flag_ensemble_docking", "w") as dock_file:
        dock_file.write("-in:file:s output_files/prepack/" + pdb[:-4] + "_prepack_0001.pdb" + "\n")
        dock_file.write("-unboundrot " + pdb + "\n\n")
        dock_file.write("#SBATCH --ntasks=" + str(cpu_dock) + "\n")
        dock_file.write("-nstruct " + str(docking_runs) + " \n\n")
        dock_file.write("-partners E_A\n")
        dock_file.write("-dock_pert 3 8\n\n")
        dock_file.write("-spin\n-detect_disulf true\n-rebuild_disulf true\n\n")
        dock_file.write("-ensemble1 protein_ensemblelist\n-ensemble2 peptide_ensemblelist\n\n")
        dock_file.write("-ex1\n-ex2aro\n\n")
        dock_file.write("-docking_low_res_score motif_dock_score\n")
        dock_file.write("-mh:path:scores_BB_BB " + rosetta_dir
                        + "/main/database/additional_protocol_data/motif_dock/xh_16_\n")
        dock_file.write("-mh:score:use_ss1 false\n")
        dock_file.write("-mh:score:use_ss2 false\n")
        dock_file.write("-mh:score:use_aa1 true\n")
        dock_file.write("-mh:score:use_aa2 true\n")
        dock_file.write("-out:path:all output_files/dock\n")
        dock_file.write("-out:suffix _ensemble_dock\n")


###############################################
# REFINE
# Returns the best pdb from the score log form docking
def check_score_dock():
    score_dic = {}
    best_pdb = ""
    score = 100000.0
    with open("output_files/dock/score_ensemble_dock.sc", "r") as score_read:
        for line in score_read:
            if not line.__contains__("SEQUENCE") and not line.__contains__("total_score"):
                score_dic[line.split()[-1]] = line.split()
    for pdb in score_dic:
        if float(score_dic[pdb][5]) < score:
            best_pdb = score_dic[pdb][-1]
            score = float(score_dic[pdb][5])
    return best_pdb


# Method: remove_dock()
# Goal: Remove files that are not highest scoring
def remove_dock(best_pdb):
    for pdb in os.listdir(os.getcwd() + "/output_files/dock/"):
        if pdb.endswith(".pdb"):
            if pdb != best_pdb + ".pdb":
                os.remove(os.getcwd() + "/output_files/dock/" + pdb)


def run_refine(pdb):
    dir_dock = rosetta_dir + "/main/source/bin/docking_protocol.mpi." + version
    make_refine_file(pdb)
    subprocess.run(["mpirun", dir_dock, "@flag_local_refine"], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)


def make_refine_file(pdb):
    with open("flag_local_refine", "w") as refine_file:
        refine_file.write("-in:file:s output_files/dock/" + pdb + "\n")
        refine_file.write("#SBATCH --ntasks=" + str(refine_runs) + "\n")
        refine_file.write("-nstruct " + str(refine_runs) + " \n\n")
        refine_file.write("-docking_local_refine\n")
        refine_file.write("-use_input_sc\n\n")
        refine_file.write("-ex1\n-ex2aro\n\n")
        refine_file.write("-out:file:fullatom\n")
        refine_file.write("-out:path:all output_files/refine\n")
        refine_file.write("-out:suffix _local_refine\n")


# Method: check_score_refine()
# Goal: Check for best scoring file with RE
def check_score_refine():
    for file in os.listdir(os.getcwd() + "/output_files/refine"):
        if file.endswith(".fasc"):
            score_file = "output_files/refine/" + file
    express = "tail -n +3 " + score_file + " | tr -s ' ' | sort -u -k6 -r | head -1"
    temp_best = subprocess.run(express, shell=True, stdout=subprocess.PIPE)
    best_pdb = temp_best.stdout.decode('utf-8').split(' ')[-1][:-1]  # Removing new line character
    return best_pdb


# Method: remove_refine()
# Goal: Remove files that are not highest scoring
def remove_refine(refine_best_pdb):
    for pdb in os.listdir(os.getcwd() + "/output_files/refine/"):
        if pdb.endswith(".pdb"):
            if pdb != refine_best_pdb + ".pdb":
                os.remove(os.getcwd() + "/output_files/refine/" + pdb)
        

###############################################


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

