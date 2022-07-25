# ACEdecoy
Repository for scripts and data produced for publication of "Designing an ACE2-derived fragment as a decoy for novel SARS-Cov-2 virus"

## Installations needed
...

## Docking Tools

##### flexauto_reosetta.py
Script to automate the running of RosettaDock v4.0. By default performs 100 relaxation runs for the spike protein. Then performs three different relax protocols for the peptide for the ensembles used for "flexible" docking. Prepacking is performed and then 5000 docking runs are performed.

To adjust number of CPU cores allocated and how many steps of relax, docking, or refinement to run, modify lines 11-22. To point to your mpi installation of Rosetta, modify line 329.

##### PDB_Tools_V3.py
Tool developed for use with PDB files (mostly TCR PDB files). Additional documentation will be posted here in the future.

##### PostCoupler.py
Tool used to re-match the number of the crystal structure reference. This is needed due to the numbering scheme required by RosettaDock.

##### Example of tools:
Align the peptide to the crystal structure of 6m0j.pdb, deleting the ACE2 atoms, and save as 6m0j_peptide1.pdb.
> python3 PDB_Tools_V3.py 6m0j_peptide1.pdb --renum2
> touch peptide.pdb
> touch protein.pdb
Copy atoms from 6m0j_peptide_renum.pdb for the peptide and protein into their respective pdb files. The following step will need to be ran on a computer with a high number of CPU cores.

Ensure you have modified flexauto_rosetta.py ask instructed above and the peptide.pdb and protein.pdb files are located in the same directory.
> python3 flexauto_rosettta.py 6m0j_peptide1.pdb

Suggest to run with nohup
> nohup python3 flexauto_rosettta.py 6m0j_peptide1.pdb &

Once the run is complete, the best pose after refinement will be copied to working directory and name printed in nohup file.

For the last step, we'll match the numbering as 6m0j for the spike protein using PostCoupler.py
> python3 PostCoupler.py 6m0j_peptide1_docked.pdb -r 6m0j.pdb -c E


## EnergyBreakdown
To use the tools, ensure that Rosetta is installed on your device and point to the directory in config.ini.

##### Rosetta_Breakdown.py
Automates the running of Rosetta's Residue Energy Breakdown script and curates into a table or heat map.

##### relax_only.py
A copy of flexauto_rosetta.py but with only relaxation of a single complex to ensure the structure is devoid of clashes and conforms to the force field used by Rosetta.

##### PDB_Tools_V3.py
Used here to renumber the structure undergoing relaxation based on Rosetta's requirements.

##### PostCoupler.py
Used here to renumber the PBD after relaxation.

##### energy_color.py
Script that can be ran in Chimera to color the interaction of the structure a csv table has been produced by Rosetta_Breakdown.py.

##### Example of tools:
This example will follow the relaxation of the PDB 6m0j.pdb.

Renumber the pdb in accordance of Rosetta.
> python3 PDB_Tools_V3.py 6m0j.pdb --renum2

Point to Rosetta install in config.ini file then before relaxation. Replace 6m0j_relaxed.pdb with file that was returned by previous step. Ensure to modify line 99 of relax_only.py to point to installation of Rosetta.
> python3 relax_only.py 6m0j_renum_relaxed.pdb

Renumber the relaxed structure back to the original numbering.
> python3 PostCoupler.py 6m0j_relax.pdb -r 6m0j.pdb -c AE

Now perform Residue Energy Breakdown. There are several ways you can analyze the structure.

Interaction energy heat map
> python3 Rosetta_Breakdown.py 6m0j_relax_renum.pdb -s E -p A -e
Interaction energy table
> python3 Rosetta_Breakdown.py 6m0j_relax_renum.pdb -s E -p A -t
Distance heat map can also be produced, for more information run...
> python3 Rosetta_Breakdown.py -h

## MDsimulation
Molecular Dynamics simulations are performed with GROMACS. To perform a run follow the commands in gromac_cmds.txt.

The tool helix_ploy.py will read in xvg output from GROMACS and gmx do_dssp and produce plots.

The folder "MDPfiles" contains the flag files for each step. These all follow the Lysozyme tutorial from GROMACS but adjusted for 100 ns runs. 

The folder "MDsimulation/Peptides" contains the plots produced from runs for paper.

## Peptides
This folder contains relevant data for each peptide and the associated Docking runs and MD simulations.
