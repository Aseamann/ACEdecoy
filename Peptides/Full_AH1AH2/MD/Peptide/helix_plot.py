import os.path
import seaborn as sns
import argparse
import re
import matplotlib.pyplot as plt
import pandas as pd
plt.style.use('seaborn-whitegrid')


def read_scount(file_in, number, multi=False):
    """
    Read in do_dssp output xvg file and organize into a dictionary for plotting with seaborn

    Parameters
    ----------
    file_in : str
        Location of DSSP output xvg file
    number : int
        Number of residues in structure
    multi : bool
        If there are multiple xvg files to collect data from

    Returns
    -------
    structure_data : dict
        Dictionary of DSSP output
    """
    files = []
    if multi:
        all_data = {}
        for file_i in os.listdir(file_in):
            if file_i.endswith(".xvg"):
                files.append(file_in + "/" + file_i)
                all_data[file_i.split(".")[0]] = {}
    else:
        files = [file_in]
    for file_i in files:
        with open(file_i, "r") as f1:
            structure_data = {}  # {structure: [0, 1]}
            for line in f1:
                if re.match("^@ s\d", line[0:4]):
                    structure_data[line.split('"')[1]] = []
                if len(structure_data.keys()) != 0 and line[0] == " ":
                    data = line.split()
                    for i in range(len(structure_data.keys())):
                        # Append number of residues in structure divided by total number of residues
                        structure_data[list(structure_data.keys())[i]].append(int(data[i + 1]) / number)
            if multi:
                all_data[file_i.split("/")[-1].split(".")[0]] = structure_data
    if multi:
        return all_data
    else:
        return structure_data


def plot_structure(data, structure, multi=False):
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
        y_data = data[structure]
        y_data = y_data[::100]
        ax = sns.lineplot(x=range(len(y_data)), y=y_data)
    plt.ylim(0, 1)
    plt.ylabel("fraction of residues in alpha-helix")
    plt.xlabel("Time (ns)")
    plt.show()


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("scount_file", help="DSSP output xvg file or folder from GROMACS do_dssp", type=str)
    parser.add_argument("-s", "--structure", help="Structure to plot", type=str)
    parser.add_argument("-n", "--number", help="Number of residues in structure", type=int)
    return parser.parse_args()


def main():
    args = parse_args()
    if os.path.isfile(args.scount_file):
        structure_data = read_scount(args.scount_file, args.number)
        plot_structure(structure_data, args.structure)
    elif os.path.isdir(args.scount_file):
        structure_data = read_scount(args.scount_file, args.number, True)
        plot_structure(structure_data, args.structure, True)


if __name__ == '__main__':
    main()
