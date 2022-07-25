import chimera
import os
import Midas
import math
from colour import Color

directory = os.getcwd()
path = directory + "6m0j_fullah1ah2_docked_orinum.pdb"
energy_table = "output.csv"
chimera.runCommand("split #0 chains")


def color_by_magnitude(value):
    white = Color("white")
    colors = list(white.range_to(Color("darkred"), 102))
    colors = colors[1:-1]
    return str(colors[int(value) - 1])


def get_mag(current, min_i, max_i):
    value = (1.0 / (max_i - min_i)) * abs(current)
    return value * 100


with open(energy_table, "r") as file:
    values = []
    ran = {}
    min_i = 10000000.0
    max_i = -10000000.0
    for line in file:
        split = line.split(",")
        if split[0] != "ACE2":
            if float(split[1]) <= 0.0:
                # 0: ACE2 AA, 1: value, 2: Spike AA
                values.append([split[0], float(split[1]), split[2][:-1]])
    for each in values:
        if each[1] > max_i:
            max_i = each[1]
        if each[1] < min_i:
            min_i = each[1]
    # Ensure most negative is shown and work towards most neg
    for each in sorted(values, key=lambda x: x[1], reverse=True):
        # print(each[0])
        # print(each[1])
        # print(each[2])
        # print(get_mag(each[1], min_i, 0))
        if each[1] <= 0.0:  # Not looking at positive values
            if int(get_mag(each[1], min_i, 0)) != 0:
                # Colors ACE2
                chimera.runCommand("color " + color_by_magnitude(get_mag(each[1], min_i, 0)) + ",r,s #0.2:"
                                   + str(each[0].split(" ")[1]))
                # Shows the side chain of ACE2
                chimera.runCommand("display #0.2:" + str(each[0].split(" ")[1]))
                # Colors Spike
                chimera.runCommand(
                        "color " + color_by_magnitude(get_mag(each[1], min_i, 0)) + ",r,s #0.1:" + str(each[2].split(" ")[1]))
                # Shows the side chain of spike
                chimera.runCommand("display #0.1:" + str(each[2].split(" ")[1]))

