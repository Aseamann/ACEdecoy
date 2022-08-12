import subprocess
import os

for pdb in os.listdir():
    if pdb.endswith(".pdb"):
        print("Running: " + pdb)
        info = ["python3", "Rosetta_Breakdown.py"]
        info.extend([pdb, "-s", "E", "-p", "A", "-t"])
        subprocess.run(info, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        os.rename("output.csv", pdb[:-4] + ".csv")
