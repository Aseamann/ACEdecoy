#!/usr/bin/python3

######################################################################
# Rosetta_Breakdown.py -- Adapted from a tool for a future pub.      #
# Copyright: Austin Seamann & Dario Ghersi                           #
# Version: 0.1                                                       #
# Last Updated: July 11th, 2022                                      #
# Goal: Perform analysis of Spike/ACE-2 iteractions with Rosetta     #
#                                                                    #
# Named arguments: -a --ab ((AB usage) Submit PDB for ab_usage       #
#                           interface scores)                        #
#                  -v --verbose ((AB usage) Verbose)                 #
#                  -s --sc ((Native) Score file produced from        #
#                           docking or refinement (or dir of .sc))   #
#                  -x --xaxis ((Native) X axis for native structure  #
#                              comparison)                           #
#                  -y --yaxis ((Native) Y axis for native structure  #
#                              comparison)                           #
#                  --ymax ((Native) Y axis value maximum)            #
#                  --ymin ((Native) Y axis value minimum)            #
#                  --xmax ((Native) X axis value maximum)            #
#                  --xmin ((Native) X axis value minimum)            #
#                  -p --heatmap_dist ((DHM) Distance Breakdown TCR   #
#                                     PDB File)                      #
#                  -d --distance ((DHM) Distance Breakdown TCR PDB   #
#                                 File)                              #
#                  -c --alpha_carbon ((DHM) Alpha carbon only)       #
#                  -e --heatmap_energy ((EB) TCR PDB File)           #
#                  -t --table ((EB) TCR PDB File)                    #
#                  -m --mhc ((DHM/EB) Changes energy breakdown to    #
#                            MHC versus peptide)                     #
#                  -f --fontsize ((DHM/EB) Adjust font size)         #
######################################################################

import argparse
import os
import sys

import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns
import subprocess
from PDB_Tools_V3 import PdbTools3
import pickle


#################
#    Global     #
#################
rosetta_dir = ""
verbose = False


#################
#    Methods    #
#################
###################
#     Heatmap     #
###################
# Global
first_aa = {}


def read_sheet(sheet_in, chain_in):
    """
    Reading energy breakdown csv file
    File: 0: Score, 1: pose_id, 2: resi1, 3: pdbid1, 4: restype1, 5: resi2, 6: pdbid2, 7: restype2, 8: fa_atr
    9: fa_rep, 10: fa_sol, 11: fa_sol_rep, 12: fa_intra_sol_xover4, 13: ik_ball_wtd, 14: fa_elec, 15: pro_close
    16: hbond_sr_bb, 17: hbond_lr_bb, 18: hbond_bb_sc, 19: hbond_sc, 20: dslf_fa13, 21: omega, 22: fa_dun
    23: p_aa_pp, 24: yhh_planarity, 25: ref, 26: rama_prepro, 27: total, 28: description

    Parameters
    ----------
    sheet_in : str
        name / location of csv file

    Returns
    -------
    chains : dict
        {'A':[['1A','GLY],['2A','SER'],...]}
    interactions : dict
        {'1A':['1A', '2A',0.0],...],'2A':[..]}
    """
    interactions = {}  # Dir of interactions in PDB, Key:
    chains = {}  # Dir of chains in file, Key: Chain_ID, Value: 0 - red_id, 1 - AA
    df = pd.read_csv(sheet_in)
    for key, value in df.iterrows():
        values = value.array
        if values[7] != "onebody":
            if values[3][-1] != values[6][-1]:  # Ensure not capturing internal interactions NEW
                interactions[values[3]] = [[values[3], values[6], values[27]]]
                if values[6] in interactions.keys():
                    interactions[values[6]].append([values[3], values[6], values[27]])
                else:
                    interactions[values[6]] = [[values[3], values[6], values[27]]]
        elif values[7] == "onebody":
            if values[3][-1] in chains.keys():
                chains[values[3][-1]].append([values[3], values[4]])
                if values[3][-1] == chain_in:
                    interactions[values[3]] = []
            else:
                first_aa[values[3][-1]] = int(values[3][:-1])
                chains[values[3][-1]] = [[values[3], values[4]]]
    # print(chains)
    # print(interactions)
    return chains, interactions


def get_peptide_inter(chains, interactions, chain_in):
    """
    Gathers each interaction and produces a dictionary with {'1A':'GLY',...,'105D':'GLY'}

    Parameters
    ----------
    chains : dict
        Dictionary produced by above method - {'A':[['1A','GLY],['2A','SER'],...]}
    interactions :dict
        Dictionary produced by above method - {'1A':['1A', '2A',0.0],...],'2A':[..]}
    chain_in : str
        Chain being compared to for interactions to the TCR

    Returns
    -------
    aa_inter : dict
        Dictionary containing the interactions to be documented in the heatmap or csv table
    """
    peptide = []
    aa_inter = {}
    for AA_Peptide in chains[chain_in]:
        peptide.append(AA_Peptide)
    for each in peptide:
        if each[0] in interactions.keys():
            if each[0] in aa_inter.keys():
                aa_inter[each[0]].append(interactions[each[0]])
            else:
                aa_inter[each[0]] = [interactions[each[0]]]
    return aa_inter


def make_peptide_table(aa_list, chain_list, csv_name, spike):
    """
    Creates the CSV table of results of the interactions

    Parameters
    ----------
    aa_list : dict
        List of interactions provided by method above
    chain_list : list
        List of chains being compared
    csv_name : str
        Name of the outputted table being produced
    spike : str
        Spike protein chain label
    """
    aa_info = {}  # Key: ChainID Value: AA as 3 letter
    for chain in chain_list:
        for each in chain_list[chain]:
            aa_info[each[0]] = each[1]
    # print(aa_list)
    with open(csv_name, "w") as t1:
        t1.write("ACE2,Interaction energy,Spike\n")
        for AA in aa_list:
            for each in aa_list[AA][0]:
                if each[0] == AA:
                    partner = each[1]
                else:
                    partner = each[0]
                if partner[-1] != spike:  # Only allows for TCR chains
                    output = aa_info[AA] + " " + str(int(AA[:-1])) + "," + str(each[2]) \
                             + "," + aa_info[partner] + " " + partner[:-1] + "\n"
                    t1.write(output)


def heatmap_info(aa_list, chain_list, chain_in):
    """
    Generates data for heatmap based on chain_in
    Return csv with dataset to create heatmap

    Parameters
    ----------
    aa_list : dict
        List of interactions provided by method above
    chain_list : list
        List of chains being compared
    chain_in : str
        Chain that's being determined what interactions it has to the TCR chains

    Returns
    -------
    file_name : str
    cdr_info : dict
    """
    aa_info = {}  # Key: ChainID Value: aa as 3 letters
    header_list = []  # Saves the header constructed: numAA
    for chain in chain_list:  # Key: '234E': 'VAL' for every AA
        for each in chain_list[chain]:
            aa_info[each[0]] = each[1]
    file_name = "temp.csv"
    with open(file_name, "w") as t2:
        for chain in chain_list:  # Chain = entire list
            if chain != chain_in:
                for aa in chain_list[chain]:  # each aa list
                    header_list.append(str(aa[0][:-1]) + chain)
                    t2.write("," + aa_info[str(aa[0][:-1]) + chain] + " " + str(aa[0][:-1]))  # Prints header ex. SER 32, VAL 50...
        t2.write('\n')
        for AA in aa_list:  # Loops though each AA in file and writes ones that partner with D or E
            t2.write(aa_info[AA] + " " + str(AA[:-1]))  # Writes column 0, peptide info
            for spike_aa in header_list:
                flag = False  # Flag for if interaction or not
                for chain_inter in aa_list[AA][0]:  # Checks each peptide_inter to tcr_aa
                    if chain_inter[0] == spike_aa:
                        t2.write("," + str(chain_inter[2]))  # Writes energy value
                        flag = True
                        break
                    elif chain_inter[1] == spike_aa:
                        t2.write("," + str(chain_inter[2]))  # Writes energy value
                        flag = True
                        break
                if not flag:
                    t2.write(",0")  # If not broken out of loop, energy 0
            t2.write("\n")
    return file_name


def heatmap(info, font_size, distance_in=0.0, distance=False, alpha_carbon=True):
    """
    Produce the heatmap with the information collected above using seaborn

    Parameters
    ----------
    info : str
        Location of CSV for heatmap creation
    chain_in : str
        Chain that's being determined what interactions it has to the spike protein
    font_size : int
        Size the user wants the font to be
    distance_in : float
        vmax
    distance : boolean
        If utilizing distance cut on plot
    """
    # Read in csv
    df = pd.read_csv(info)
    if not distance:
        # Remove rows with only zero values
        df = df[(df.sum(axis=1) != 0)]
        # Remove columns that don't contain data in-between them
        df = df.loc[:, (df != 0).any(axis=0)]
    if distance:
        flag = False
        for i in df.columns:
            if flag:
                column_list = df[i].values.tolist()
                column_list = [i for i in column_list if i != 10000]
                if len(column_list) >= 1:
                    pass
                else:
                    df = df.drop(i, axis=1)
            flag = True
        flag = False
        print(df)
        for i, row in df.iterrows():
            if flag:
                row_list = list(row)
                row_list = [i for i in row_list[1:] if i != 10000]
                if len(row_list) >= 1:
                    pass
                else:
                    df = df.drop(i)
            flag = True
    print(df)
    # Read in x and y axis labels
    y_axis_labels = list(df.iloc[:,0])
    x_axis_labels = list(df.iloc[0:, :])[1:]
    # Drop first column
    df = df.iloc[:, 1:]
    # Annotation labels
    if not distance:
        labels = df.iloc[:, :].applymap(lambda v: str(v) if v != 0 else '')
        ax = sns.heatmap(df, vmax=0, linewidths=.2, linecolor="grey", xticklabels=x_axis_labels,
                         yticklabels=y_axis_labels,
                         cmap=sns.cubehelix_palette(start=2, rot=0, reverse=True, dark=0, light=1, as_cmap=True),
                         annot=labels, annot_kws={"fontsize": font_size}, fmt='')
        ax.set_yticklabels(y_axis_labels, size=font_size)
        ax.set_xticklabels(x_axis_labels, size=font_size)
        plt.title("Pairwise Rosetta Energy Breakdown")
    if distance:
        labels = df.iloc[:, :].applymap(lambda v: str(f'{v:.3f}') if v < 10000 else '')
        ax = sns.heatmap(df, vmax=distance_in, vmin=0, linewidths=.2, linecolor="grey", xticklabels=x_axis_labels,
                         yticklabels=y_axis_labels,
                         cmap=(sns.cubehelix_palette(start=2, rot=0, reverse=True, dark=0, light=1)),
                         annot=labels, annot_kws={"fontsize": font_size}, fmt='')
        ax.set_yticklabels(y_axis_labels, size=font_size)
        ax.set_xticklabels(x_axis_labels, size=font_size)
        plt.title("Pairwise Distance (Angstroms)")
        if alpha_carbon:
            plt.title("Alpha Carbon Pairwise Distance (Angstroms)")
    plt.xlabel("Spike Protein", size=font_size)
    y_label = "ACE2 Protein/Peptide"
    plt.ylabel(y_label, size=font_size)
    plt.show()


def interface_heatmap(interface_breakdown, partner, font_size):
    """
    Controls the helper method to produce heatmap of residue energy breakdown results

    Parameters
    ----------
    interface_breakdown : str
        Output of interface breakdown
    partner : str
        Chain id of partner protein/peptide
    font_size : int
        Size the user wants the font to be
    """
    chain_list, inter_list = read_sheet(interface_breakdown, partner)
    aa_inter = get_peptide_inter(chain_list, inter_list, partner)
    info = heatmap_info(aa_inter, chain_list, partner)
    heatmap(info, font_size)
    os.remove(info)


def get_contacts(pdb, distance, alpha_carbon, spike, partner):
    """
    Calculate distances between spike and partner - save atoms within cutoff distance

    Parameters
    ----------
    pdb : str
        Location of PDB file
    distance : float
        Cutoff of collected values
    alpha_carbon : boolean
        Distance from alpha_carbon - else distance from closest atom in aa
    chain_in : chain_in
        partnering chain

    Returns
    -------
    contacts : dict
        {AA comp_num: {partner AA comp_num: distance}}
    cdr_info : dict
        Dictionary containing information about CDR regions of each chains
    cdr_aa : dict
        Dictionary containing information about CDR regions for each chain
    all_aa : dict
        Store list of residues and resi num for each amino acid in each chain
    """
    tool = PdbTools3(pdb)
    spike_atoms = tool.get_atoms_on_chain(spike)  # Collect spike chain atoms
    partner_atoms = tool.get_atoms_on_chain(partner)  # Collect partner chain atoms
    all_aa = {spike: {}, partner: {}}  # Store list of residues and resi num for each amino acid in each chain
    # all_atoms = {spike: spike_atoms, partner: partner_atoms}
    # Loop through all atoms
    for chain in [spike_atoms, partner_atoms]:
        for atom in chain:
            if atom["comp_num"] not in all_aa[atom["chain_id"]].keys():  # Store all_aa info {chain: {202: THR}}
                all_aa[atom["chain_id"]][atom["comp_num"]] = atom["atom_comp_id"]
    if alpha_carbon:
        spike_atoms = [atom for atom in spike_atoms if atom["atom_id"] == "CA"]
        partner_atoms = [atom for atom in partner_atoms if atom["atom_id"] == "CA"]
    contacts = {spike: {}}  # chain_id: {AA comp_num: {partner AA comp_num: distance}}
    flag = False  # Once we've looped once
    for atom_1 in spike_atoms:
        for atom_2 in partner_atoms:
            if not flag:
                if atom_2["comp_num"] not in all_aa[atom_2["chain_id"]]:
                    all_aa[atom_2["chain_id"]][atom_2["comp_num"]] = atom_2["atom_comp_id"]
            # Calculate distance between each partner
            euc_dist = tool.euclidean_of_atoms(atom_1["atom_num"], atom_2["atom_num"])
            if euc_dist <= distance:  # if within threshold
                if atom_1['comp_num'] in contacts[spike].keys():  # if cdr aa already documented
                    if atom_2['comp_num'] in contacts[spike][atom_1['comp_num']]:  # If cdr aa partner found
                        if euc_dist < contacts[spike][atom_1['comp_num']][atom_2['comp_num']]:
                            contacts[spike][atom_1['comp_num']][atom_2['comp_num']] = euc_dist
                    else:
                        contacts[spike][atom_1['comp_num']][atom_2['comp_num']] = euc_dist
                else:  # if not previous aa to aa contact documented
                    contacts[spike][atom_1['comp_num']] = {atom_2['comp_num']: euc_dist}
        flag = True
    return contacts, all_aa  # Chain_id: {aa comp_num: {partner AA comp_num: distance}}


def distance_heatmap_info(contacts, all_aa, spike, partner):
    """
    Produce the heatmap csv information table

    Parameters
    ----------
    contacts : dict
        {AA comp_num: {partner AA comp_num: distance}}
    all_aa : dict
        Store list of residues and resi num for each amino acid in each chain
    spike : str
        Chain id of spike protein
    partner : str
        Chain id of partner protein

    Returns
    -------
    file_name : str
        Name of csv produced
    """
    file_name = "temp.csv"
    with open(file_name, "w") as f1:
        header_list = []  # Keeps tracks of position of aa in header list ex. 101E
        for num in all_aa[spike]:  # Loop over each position in cdr
            # Write Three Letter + Resi Num
            if num in all_aa[spike].keys():
                header_list.append(str(num) + spike)
                f1.write("," + all_aa[spike][num] + " " + str(num))
        f1.write("\n")
        for num in all_aa[partner]:
            f1.write(all_aa[partner][num] + " " + str(num))
            for aa in header_list:
                if int(aa[:-1]) in contacts[aa[-1]].keys():
                    if num in contacts[aa[-1]][int(aa[:-1])].keys():  # Contacts[chain][num]
                        f1.write("," + str(contacts[aa[-1]][int(aa[:-1])][num]))
                    else:
                        f1.write("," + str(10000))
                else:
                    f1.write("," + str(10000))
            f1.write("\n")
    return file_name


def interface_heatmap_dist(pdb, distance, alpha_carbon, spike, partner, font_size):
    """
    Controls the methods used for interface heatmap with distance cutoff

    Parameters
    ----------
    pdb : str
        Location of PDB file
    distance : float
        Distance cutoff
    alpha_carbon : boolean
        If using alpha_carbon distance
    partner : str
        Chain id of partner protein/peptide
    font_size : int
    """
    pickle_name = pdb[:-1] + "_contacts.pickle"
    pickle2_name = pdb[:-1] + "_all_aa.pickle"
    if pickle_name not in os.listdir() and pickle2_name not in os.listdir():
        contacts, all_aa = get_contacts(pdb, distance, alpha_carbon, spike, partner)  # Collect contacts
        with open(pickle_name, "wb") as f:
            pickle.dump(contacts, f)
        with open(pickle2_name, "wb") as f1:
            pickle.dump(all_aa, f1)
    else:
        with open(pickle_name, 'rb') as f:
            contacts = pickle.load(f)
        with open(pickle2_name, 'rb') as f1:
            all_aa = pickle.load(f1)
    temp_file = distance_heatmap_info(contacts, all_aa, spike, partner)  # Create csv for heatmap
    heatmap(temp_file, font_size, distance, True, alpha_carbon)  # Create heatmap
    # os.remove(temp_file)


def peptide_table(energy_breakdown, partner):
    """
    Controller method to generate energy breakdown table

    Parameters
    ----------
    energy_breakdown : str
        Location fo energry breakdown output
    partner : str
        Chain id of partnering peptide/protein
    """
    chain_list, inter_list = read_sheet(energy_breakdown, partner)
    aa_inter = get_peptide_inter(chain_list, inter_list, partner)
    make_peptide_table(aa_inter, chain_list, "output.csv", partner)


def tsv_to_csv(tsv_in):
    """
    Convert tsv to csv

    Parameters
    ----------
    tsv_in : str
        Location of tsv file

    Returns
    -------
    new_name : str
        Updated location of file, now csv
    """
    name_in = tsv_in
    # Rename to same name but replace file extension to .csv
    new_name = "".join(tsv_in.split(".")[:-1]) + ".csv"
    submit_name = name_in.replace(" ", "\\ ")
    # Convert .out (tsv) to .csv
    cmd = "tr -s ' ' < " + submit_name + " | tr ' ' ','"
    new_file = subprocess.run([cmd], shell=True, stdout=subprocess.PIPE)
    os.rename(name_in, new_name)
    with open(new_name, "w") as f1:
        f1.write(new_file.stdout.decode('utf-8')[1:])
    return new_name


def run_breakdown(pdb_in):
    """
    Run Rosetta Residue Energy Breakdown program if PDB submitted for heatmap or energy breakdown table

    Parameters
    ----------
    pdb_in : str
        Location of PDB file

    Returns
    -------
    convert_out : str
        Location of the outputted csv results
    """
    # Update program location information
    global rosetta_dir
    # In and out files
    in_file = "-in:file:s " + pdb_in
    file_name = "energy_breakdown.out"
    out_file = "-out:file:silent " + file_name
    # Run process
    subprocess.run([rosetta_dir, in_file, out_file], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    # Convert tsv to csv
    file_name = os.getcwd() + "/" + file_name
    convert_out = tsv_to_csv(file_name)
    return convert_out


def rosetta_binary(program_in):
    """
    Determine what binaries are built for rosetta based on each program

    Parameters
    ----------
    program_in : str
        Location of Rosetta
    """
    # Update program location information
    global rosetta_dir
    with open("config.ini", "r") as f1:
        for line in f1:
            if line[:11] == "rosetta_loc":
                rosetta_dir = line[:-1].split("=")[1][1:-1]
    if not rosetta_dir.endswith("/"):
        rosetta_dir += "/"
    programs = []
    for program in os.listdir(rosetta_dir + "main/source/bin"):
        if program.startswith(program_in):
            programs.append(program)
    for program in sorted(programs, key=len):
        if program.startswith(program_in + ".mpi"):
            rosetta_dir += "main/source/bin/" + program
            break
        elif not program.startswith(program_in + ".default"):
            rosetta_dir += "main/source/bin/" + program


####################
#     Controls     #
####################
def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("pdb", help="Location of pdb", type=str)
    parser.add_argument("-s", "--spike", help="Label for spike protein", type=str)
    parser.add_argument("-p", "--partner", help="Label for partnering protein/peptide", type=str)
    parser.add_argument("-i", "--heatmap_dist", help="(DHM) Produce distance heatmap", default=False,
                        action="store_true")
    parser.add_argument("-d", "--distance", help="(DHM) Distance for heatmap | Default 4.5", type=float, default=4.5)
    parser.add_argument("-c", "--alpha_carbon", help="(DHM) Alpha Carbon Only", action="store_true", default=False)
    parser.add_argument("-e", "--heatmap_energy", help="(EB) Produce energy breakdown heatmap", default=False,
                        action="store_true")
    parser.add_argument("-t", "--table", help="(EB) Produce energy breakdown table", default=False, action="store_true")
    parser.add_argument("-f", "--fontsize", help="(DHM/EB) Adjust font size | Default 12", type=int, default=12)
    return parser.parse_args()


def main():
    args = parse_args()
    if args.heatmap_energy or args.table:  # Generate heatmap
        if args.pdb.endswith(".pdb"):  # Run energy breakdown if PDB submitted
            rosetta_binary("residue_energy_breakdown")  # Set binary
            breakdown_file = run_breakdown(args.pdb)
            if args.heatmap_energy:  # Updated
                interface_heatmap(breakdown_file, args.partner, args.fontsize)
            if args.table:
                peptide_table(breakdown_file, args.partner)
            os.remove(breakdown_file)
    # Heatmap for distance
    if args.heatmap_dist:
        if args.pdb.endswith(".pdb"):
            interface_heatmap_dist(args.pdb, args.distance, args.alpha_carbon, args.spike, args.partner, args.fontsize)
        else:
            print("Provide PDB File")


if __name__ == '__main__':
    main()
