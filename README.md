**Structural replacement for the design of protein - small molecule binding** 

two options for modeling:
* dock to sequence 
* glycine shaved docking 


**Required Software:**

* Python 3.8+
* A local installation of [PyRosetta](https://www.pyrosetta.org/home)
* [numPy](https://numpy.org/)
* [pandas](https://pandas.pydata.org/)


**General Usage:**

We recommend generating small molecule conformers using [SM_ConfGen](https://github.com/ajfriedman22/SM_ConfGen/tree/main).  

This repository builds on the Structural Replacement method developed by Jordan Wells as part of a 2021 Rosetta summer internship. Many thanks to Jordan for his contributions. His code can be found here and is an excellent reference: https://github.com/jordantwells42/structural-replacement. The portions of this method that directly use Jordan's scripts have his adapted his instructions for use below.

Each script is ran with the path to the config file as an argument, such as

`python path/to/a_script.py path/to/conf.txt`

**Config File:**

The config file contains all of the information needed for the scripts to run.

Information in the [DEFAULT] category is used by all scripts, and each script has its own respective category.

*REQUIRED* options need to be filled in to fit your own specific project’s needs, while *OPTIONAL* options do NOT need to be filled in and can be left as is. The default settings in the optional options are what I have found to work the best.

Additionally, a lot of the options in the config files are paths to locations on your computer. Importantly these paths are in relation to where your scripts are located, NOT the config file.

**Running the Scripts:**

Once the config file has been set up with all of the necessary information. The necessary commands and a brief outline of what each script does is outlined below. If you would like more information on what each script does or how to use them, that can be found [here](https://docs.google.com/document/d/1NEq-mbIoxclpstKW4C55wvxyhdNPPFbA7jrmUidYLdk/edit?tab=t.0).

`python create_table.py conf.txt`
* This will create a spreadsheet from the ligand sdf files provided in the *MoleculeSDFs* option, create Rosetta .params files to make them Rosetta-readable, and create .pdb files for each ligand
* If you want to have several conformers for each ligand, simply have all of the conformers for each ligand in one file and pass that to MoleculeSDFs

*Manual input of atom alignments*
* Here you will fill in the resulting “Molecule Atoms” and “Target Atoms” columns in the generated spreadsheet by adding the atom labels for each that can be found with a software such as PyMOL by clicking on the atoms.
* For example, if I wanted to align indole to Tryptophan, I’d go into PyMOL and find a substructure they share in common (such as the six-membered ring), choose corresponding atoms on that substructure, and list the atom labels .
* This would ultimately look like a “C1-C5-C7” in the “Molecule Atoms” columns and a “"CD2-CZ2-CZ3” in the “Target Atoms” column (by chooosing three correpsonding atoms on the six-membered ring)
* Copy the row of alignment information as many times as you would like to sample alignments in the next step

*Grading conformer alignments*
* This step samples molecule alignments to identify possible alignments that do not clash with the protein structure. There are two different scripts you can run depending on which method you need
* For both methods, in config file set IncludeSC = True

`python grade_conformers_docked_to_sequence.py conf.txt` 
This method allows you to find molecule conformers and alignments that are compatible with the protein backbone and side chains of a set amino acid sequence, usually one that has been identified as binding your molecule in a low-affinity screen. 

`python grade_conformers_glycine_shaved_docking.py conf.txt` 
This method allows you to dock the ligand conformers into a glycine-shaved ligand binding pocket, essentially considering only the protein backbone at defined positions. This is useful if you want to map where your ligand could potentially fit in the pocket. 


`python rosetta_design_score_passing_pdbs.py conf.txt`
* Before running this script you will need to update your config file with the numbers of the pdb files you want to design around, found in the folder "pass_score_repacked". 





Protocol Capture 

PyRosetta-4 2021 [Rosetta PyRosetta4.Release.python38.mac 2021.26+release.b308454c455dd04f6824cc8b23e54bbb9be2cdd7 2021-07-02T13:01:54] retrieved from: http://www.pyrosetta.org
Python 3.8.13
numpy 1.21.5
pandas 1.4.2


