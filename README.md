**Structural Replacement Biosensor Design** 

Leonard, A.C. et al. Computational design of dynamic biosensors for emerging synthetic opioids. BioRxiv. https://www.biorxiv.org/content/10.1101/2025.05.15.654300v1.full.pdf

**Required Software:**

* Python 3.8+
* A local installation of [PyRosetta](https://www.pyrosetta.org/home)
* [numpy](https://numpy.org/)
* [pandas](https://pandas.pydata.org/)


**General Usage:**

The Structural Replacement computational docking protocol models protein binding to a small molecule ligand by superimposing the novel ligand on an existing solved protein-ligand crystal structure. By assuming that the novel ligand will dock in the same binding pocket as the previous ligand and using a similar mechanism, we dramatically reduce the required search space during ligand docking, providing a starting point from which to sample plausible ligand alignments during molecular design. 

Required files: 
1. A config.txt file to control docking and design settings 
2. PrePDBFile: The solved protein-ligand structure to use as a template for docking. This PDB file should be cleaned of all non-structural information and anisotropic temperature factor (ANISOU) lines. Any water molecules essential for ligand alignment can be retained. The ligand to be used as the template for alignment should be positioned as residue 1 on its own chain X. 
3. PostPDBFile: The same file as PrePDFFile, but edited to remove the template ligand. 
4. A params file for any ligands in the PDB files. This will usually be the small molecule ligand that you are replacing during docking. This is required because ligands are not Rosetta-readable without supplying this documentation. Instructions for generating params files and saving them in your PyRosetta database can be found here (https://www.pyrosetta.org/documentation/obtaining-and-preparing-ligand-pdb-files)
5. Generated ligand conformers of your desired ligand to bind. We recommend generating small molecule conformers by TREMD using [SM_ConfGen](https://github.com/ajfriedman22/SM_ConfGen/tree/main), saving the conformers in a single .sdf file. 
6. A Rosetta-relaxed version of PostPDBFile and a resfile if performing Rosetta design. 


This method has two options for modeling protein-ligand docking: 
* Dock to sequence: use this option when you have a specific protein sequence that you want to dock to, such as a low-affinity binding sequence identified experimentally that you wish to optimize. This method allows you to chose positions in the binding pocket to mutate from the wild-type sequence to match your known binding sequence. 

* Glycine shaved docking: use this option when you have no sequence information to guide docking. This method allows you to select positions in the binding pocket to shave to glycine, creating space to sample all plausible ligand alignments not constrained by the protein backbone. 


**Instructions to run:** 

Each script is ran with the path to the config file as an argument, such as

`python path/to/a_script.py path/to/conf.txt`

Config File: the config file contains all of the information needed for the scripts to run.

Information in the [DEFAULT] category is used by all scripts, and each script has its own respective category. Ensure that all file paths are in relation to where your scripts are located, NOT the config file.

*REQUIRED* options need to be filled in to fit your own specific project’s needs, while *OPTIONAL* options do NOT need to be filled in and can be left as is. The default settings in the optional options are what we have found to work well. 

*DOCKED TO SEQUENCE* and *GLYCINE SHAVED DOCKING* are both specified in the [grade_conformers] section. Choose the method you wish to run and input the information required. Leave the section blank if you are not running that method. 

**Scripts to run:**
Run each script separately and verify the output is as intended before proceeding to the next step. 

`python path/to/create_table.py path/to/config.txt`
* This will create a the spreadsheet Ligands.csv from the ligand sdf files provided in the *MoleculeSDFs* option, create Rosetta .params files to make them Rosetta-readable, and create .pdb files for each ligand (note: you can rename this spreadsheet in the config fig if you choose)
* If you want to have several conformers for each ligand, simply have all of the conformers for each ligand in one file and pass that to MoleculeSDFs
* You can align the same conformer set in different ways by listing the same comformer set multiple times in the MoleculeSDFs line of the config file create table section. See the file config_example.txt in the folder files_for_PYR1_docking as an example. 

*Manual input of atom alignments*
* Here you will fill in the resulting “Molecule Atoms” and “Target Atoms” columns in the generated spreadsheet Ligands.csv by adding the atom labels for each. Atom labels can be found using PyMOL by clicking on the atoms.
* For example, if I wanted to align indole to Tryptophan, I’d go into PyMOL and find a substructure they share in common (such as the six-membered ring), choose corresponding atoms on that substructure, and list the atom labels .
* This would ultimately look like a “C1-C5-C7” in the “Molecule Atoms” columns and a “"CD2-CZ2-CZ3” in the “Target Atoms” column (by chooosing three correpsonding atoms on the six-membered ring)
* Copy the row of alignment information as many times as you would like to sample alignments in the next step. When performing the dock to sequence method, you may need to sample several hundred total alignments to find one that matches your input sequence. 

*Grading conformer alignments*
* This step samples molecule alignments to identify possible alignments that do not clash with the protein structure. Choose the script to run depending on the method you are using. 

`python grade_conformers_docked_to_sequence.py conf.txt` 
* This method allows you to find molecule conformers and alignments that are compatible with the protein backbone and side chains of a set amino acid sequence, usually one that has been identified as binding your molecule in a low-affinity screen. This step generates PDB files containing protein-ligand structures binding the your input amino acid sequence. 

`python grade_conformers_glycine_shaved_docking.py conf.txt` 
* This method allows you to dock the ligand conformers into a glycine-shaved ligand binding pocket, essentially considering only the protein backbone at defined positions. This is useful if you want to map where your ligand could potentially fit in the pocket. 

`python rosetta_design_score_passing_pdbs.py conf.txt`
* This script allows you to perform Rosetta sequence design on any passing PDBs you select. Use this script instead of running design on individual PDB files to easily handle ligand molecules without needing to further define them in Rosetta. 
* Update the [rosetta_design] section of the Config file to include a Rosetta-relaxed version of PostPDBFile, a Rosetta resfile specifing the positions to mutation and the allowable mutations, and the PDB file numbers of the passing alignments you want to design around. The passing PDB files can be found in the folder "pass_score_repacked" and can be viewed in PyMOL. 

If you have Rosetta installed, you can relax a pdb file with the following command, choosing the correct relax release for your computer operating system. 
`path/to/Rosetta/main/source/bin/relax.macosclangrelease -database ~path/to/Rosetta/main/database -s filename.pdb -out:suffix _relaxed @relax_flags.opts`




Protocol Capture 

PyRosetta-4 2021 [Rosetta PyRosetta4.Release.python38.mac 2021.26+release.b308454c455dd04f6824cc8b23e54bbb9be2cdd7 2021-07-02T13:01:54] retrieved from: http://www.pyrosetta.org \
Python 3.8.13 \
numpy 1.21.5 \
pandas 1.4.2


This repository "Structural Replacement Biosensor Design" was written by Alison Leonard and Jordan Wells. This repository builds on the Structural Replacement method developed by Jordan Wells as part of a 2021 Rosetta summer internship. His code can be found here and is an excellent reference: https://github.com/jordantwells42/structural-replacement.

