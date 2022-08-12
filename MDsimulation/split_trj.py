import argparse
import os
import shutil
import subprocess


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("pdb", help="trj.pdb file with several poses of md run", type=str)
    parser.add_argument("-f", "--frames", help="distribution of frames to pull", type=int)
    return parser.parse_args()


def main():
    args = parse_args()
    current_pose = ""
    pose_num = 0
    for line in reversed(open(args.pdb).readlines()):
        if "MODEL" in line:
            total_pose = int(line[-5:])
            break
    with open(args.pdb, "r") as f1:
        for line in f1:
            current_pose += line  # Write every line till end of model
            if line.find("ENDMDL") != -1:
                pose_num += 1
                if pose_num % args.frames == 0 or pose_num == 1:
                    pose_name = args.pdb[:-4] + "_" + str(pose_num).zfill(len(str(total_pose))) + ".pdb"
                    with open(pose_name, "w") as f2:
                        f2.write(current_pose)
                current_pose = ""


if __name__ == "__main__":
    main()
