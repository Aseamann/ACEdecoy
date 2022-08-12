import argparse
import subprocess
import sys
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from PDB_Tools_V3 import PdbTools3
from shutil import copyfile
plt.style.use('seaborn-whitegrid')


def process_pdb(pdb, structure, rmsd):
    """
    Splits pdb into overwriting pdbs and submits command to have dssp
    process secondary structure of each pose.

    Parameters
    ----------
    pdb : str
        Location of structural information in pdb format
    structure : str
        Secondary structure to process
    rmsd : boolean
        Calculate superimposed alpha carbon rmsd and return under structure_data

    Returns
    -------
    output : list
        List of secondary structure count per pose (or rmsd if calculating rmsd)
    """
    tool = PdbTools3()
    current_pose = ""  # Going to contain string for each pose of pdb
    tmp_pdb_name = "trj_pose_tmp.pdb"
    reference_structure = "trj_pose_reference.pdb"
    # Collection Information
    structure_data = []  # Number of residues in secondary structure per pose
    total_res = 0  # Total number of residues in structure
    res_ss = 0  # Residues in secondary structure
    rmsd_flag = False
    with open(pdb, "r") as f1:
        for line in f1:
            current_pose += line  # Write every line till end of model
            # Process pose
            if line.find("ENDMDL") != -1:
                # Write pose to tmp pdb
                with open(tmp_pdb_name, "w") as f2:
                    f2.write(current_pose)
                    current_pose = ""  # Clear previous pose
                # If running secondary structure classification
                if not rmsd:
                    # Submit to mkdssp
                    express = "mkdssp -i " + tmp_pdb_name
                    dssp_output = subprocess.run([express], shell=True, stdout=subprocess.PIPE)
                    dssp_output = dssp_output.stdout.decode('utf-8')
                    flag = False
                    for line_ in dssp_output.split("\n")[:-1]:
                        if not flag and line_.find("#") != -1:
                            flag = True
                        if flag:
                            total_res += 1
                            if line_[16] == structure:
                                res_ss += 1
                    # Write to save data (add additional residue due to dssp skipping c-term resi)
                    structure_data.append(res_ss/(total_res + 1))
                    # Reset
                    res_ss = 0
                    total_res = 0
                # If calculating rmsd
                else:
                    if not rmsd_flag:  # Save first pose for calculating rmsd against
                        copyfile(tmp_pdb_name, reference_structure)
                        rmsd_flag = True
                    else:
                        tool.set_file_name(reference_structure)
                        chains = tool.get_chains()
                        # Collect RMSD
                        tool.superimpose(tmp_pdb_name, chains, chains, "tmp.pdb")
                        tool.set_file_name(tmp_pdb_name)
                        structure_data.append(tool.rmsd("tmp.pdb", chains, chains, ca=True, mute=True))
    return structure_data


def plot_structure(data, structure, multi=False):
    """
    Plot secondary structure as percent of residues in secondary structure

    Parameters
    ----------
    data
    structure
    multi

    Returns
    -------

    """
    fig, ax = plt.subplots()
    if multi:
        data_dict = {}
        for file_i in data:
            y_data = data[file_i][structure]
            y_data = y_data[::100]
            data_dict[file_i] = y_data
        data_df = pd.DataFrame(data_dict)
        sns.lineplot(data=data_df)
        # plt.legend(labels=[list(data.keys())])
    else:
        y_data = data
        y_data = y_data[::10]
        ax = sns.lineplot(x=range(len(y_data)), y=y_data)
    plt.ylim(0, 1)
    structure_label = {"H": "alpha helix", "B": "beta bridge", "E": "strand",
                       "G": "helix-3", "I": "helix-5", "T": "turn", "S": "bend"}
    plt.ylabel("fraction of residues in " + structure_label[structure])
    plt.xlabel("Time (ns)")
    # plt.show()
    plt.savefig("dssp_results.png")


def plot_rmsd(data, multi=False):
    fig, ax = plt.subplots()
    if multi:
        sys.exit()
    else:
        y_data = data
        y_data = y_data[::10]
        ax = sns.lineplot(x=range(len(y_data)), y=y_data)
    # plt.ylim(0,1)
    plt.xlabel("Time (ns)")
    plt.ylabel("alpha carbon rmsd")
    plt.savefig("rmsd_results.png")


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("pdb", help="MD output PDB containing several models in single pdb", type=str,
                        default="trj.pdb")
    parser.add_argument("-s", "--structure", help="Secondary structure to plot (ex. H, B, ...)", type=str, default="H")
    parser.add_argument("-r", "--rmsd", help="Calculate superimposed alpha carbon rmsd of each pose",
                        action="store_true", default=False)
    return parser.parse_args()


def main():
    args = parse_args()
    # Collect data to plot
    if args.rmsd:
        rmsd_data = process_pdb(args.pdb, "-", args.rmsd)
        plot_rmsd(rmsd_data)
    else:
        structure_data = process_pdb(args.pdb, args.structure, args.rmsd)
        # Plot data
        plot_structure(structure_data, args.structure)


if __name__ == "__main__":
    main()
