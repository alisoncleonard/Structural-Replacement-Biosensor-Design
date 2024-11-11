import sys
import os
import time

import pandas as pd

from configparser import ConfigParser
import argparse



def csv_to_df_pkl(csv_file_name, pkl_file_name, auto, path_to_conformers, pose, target_res, lig_res_num):
    """
    Automatically converts the CSV made by create_table into a DataFrame with all of the necessary
    information

    Creates a lot of helper functions to do so as required by Pandas

    """
    df = pd.read_csv(f"{csv_file_name}", index_col = False)
    if pkl_file_name == None:
        pkl_file_name = f"{csv_file_name[:-4]}.pkl"

    # Creating the "file stem", the location of where the ligands are
    def create_file_stem(name):
        return f"{name}/{name}"
    
    df["Molecule File Stem"] = df["Molecule ID"].apply(create_file_stem)


    # Turning the conformer range into a tuple
    def split_crange(name):
        lst = name.split("_")
        return (lst[0], lst[1])

    df["Conformer Range"] = df["Conformer Range"].apply(split_crange)


    # Auto generating alignments if specified to do so
    if auto:
        print("Auto Generating Alignments")
        for i, row in df.iterrows():
            print(i)
            print(f"{i+1}/{len(df)}", end = " ")
            lig = Pose()
            mol_id = row["Molecule ID"]
            conf_num = 1
            res_set = pyrosetta.generate_nonstandard_residue_set(lig, params_list = [f"{path_to_conformers}/{mol_id}/{mol_id}.params"])
            pyrosetta.pose_from_file(lig, res_set, f"{path_to_conformers}/{mol_id}/{mol_id}_{conf_num:04}.pdb")

            molecule_atoms, target_atoms = alignment.auto_align_residue_to_residue(lig, lig.residue(lig_res_num), target_res)
            df.loc[i, "Molecule Atoms"] = "-".join(molecule_atoms)
            df.loc[i, "Target Atoms"] = "-".join(target_atoms)


    # Turning molecule and target atoms into tuples
    def split_alabels(name):
        if name == "default":
            return ("CD2", "CZ2", "CZ3")

        lst = name.split("-")
        return tuple([str(e) for e in lst])

    df["Molecule Atoms"] = df["Molecule Atoms"].apply(split_alabels) 
    df["Target Atoms"] = df["Target Atoms"].apply(split_alabels)


    # Serializing the resultant DataFrame
    df.to_pickle(pkl_file_name, protocol = 4)

def align_to_residue_and_check_collision(pose, res, path_to_conformers, df, pkl_file,
         bin_width = 1, vdw_modifier = 0.7, include_sc = False, lig_res_num = 1):
    """
    Aligns and then checks for collisions

    Arguments
    poses: input poses of interest
    reses: input residues to map onto
    path_to_conformers: where conformers are stored
    pkl_files: pkl_files that contain all the necessary info as generated in conformer_prep
    bin_width: grid width of the collision grid
    vdw_modifier: by what factor to multiply pauling vdw radii by in the grid calculation
    include_sc: whether to do just backbone or include sc atoms

    Writes to the provided pkl files with conformers that are accepted in a column called
    "Accepted Conformers"
    """

    pmm = pyrosetta.PyMOLMover()
    pmm.keep_history(True)

    # Precalculating the protein collision grid and storing the results in a hashmap (grid.grid)
    grid = collision_check.CollisionGrid(pose, bin_width = bin_width, vdw_modifier = vdw_modifier, include_sc = include_sc)


    t0 = time.time()

    all_accepted = []
    accepted_conformations = []
    every_other = 0
    total_confs = 0


    for pose_info in conformer_prep.yield_ligand_poses(df = df, path_to_conformers = path_to_conformers, post_accepted_conformers = False, ligand_residue = lig_res_num):
        
        if not pose_info:
            print(f"{conf.name}, {conf.id}: {len(accepted_conformations)/conf.conf_num}")
            all_accepted.append(accepted_conformations)
            accepted_conformations = []
            continue

        # Grab the conformer from the generator
        conf = pose_info
        
        # Perform alignment
        conf.align_to_target(res)

        # Check for collision
        does_collide = conf.check_collision(grid)

        total_confs += 1

        if not does_collide:
            accepted_conformations.append(conf.conf_num)
            
            if every_other % 25 == 0:
                conf.pose.pdb_info().name(f"{conf.name}, {conf.id}")
                pmm.apply(conf.pose)

            every_other += 1


    print(f"\n\n---Output, {pkl_file}---")
    #print(f"List of Acceptances: {all_accepted}")
    print(f"\nNumber of Ligands: {len(all_accepted)}")
    ligands_accepted = len([e for e in all_accepted if e])
    print(f"Number of Ligands Accepted: {ligands_accepted}")
    print(f"Proportion of Ligands Accepted: {ligands_accepted/len(all_accepted)}")

    total_accepted = sum([len(e) for e in all_accepted])
    print(f"\nNumber of Conformers: {total_confs}")
    print(f"Number of Conformers Accepted: {total_accepted}")
    print(f"Proportion of Conformers Accepted: {total_accepted/total_confs}")

    tf = time.time()
    print(f"\nTime taken: {(tf - t0)/60} minutes")
    print(f"Conformers per minute: {total_confs/(tf-t0)*60}")
    
    df["Accepted Conformers"] = all_accepted
    df.to_pickle(pkl_file, protocol = 4)


def main(argv):
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("config_file",
                        help = "your config file",
                        default = "my_conf.txt")
    
    if len(argv) == 0:
        print(parser.print_help())
        return

    args = parser.parse_args(argv)

    # Parsing config file
    config = ConfigParser()
    config.read(args.config_file)
    default = config["DEFAULT"]
    spec = config["grade_conformers"]
    sys.path.append(default["PathToPyRosetta"])
    auto = default["AutoGenerateAlignment"] == "True"
    
    # Importing necessary dependencies
    global pyrosetta, Pose, alignment, conformer_prep, collision_check

    import pyrosetta
    from pyrosetta.rosetta.core.pose import Pose
    import alignment
    import conformer_prep
    import collision_check

    pyrosetta.init("-mute all")  

    params_list = default["ParamsList"]

    # Reading in Pre and Post PDBs
    print("Reading in Pre and Post PDBs")
    if len(params_list) > 0:
        pre_pose = Pose()
        res_set = pyrosetta.generate_nonstandard_residue_set(pre_pose, params_list = params_list.split(" "))
        pyrosetta.pose_from_file(pre_pose, res_set, default["PrePDBFileName"])

        post_pose = Pose()
        res_set = pyrosetta.generate_nonstandard_residue_set(post_pose, params_list = params_list.split(" "))
        pyrosetta.pose_from_file(post_pose, res_set, default["PostPDBFileName"])
    else:
        pre_pose = pyrosetta.pose_from_pdb(default["PrePDBFileName"])
        post_pose = pyrosetta.pose_from_pdb(default["PostPDBFileName"])

    
    # ACL insert code for manual glycine shave of everything we want to only look at backbone for 

    # make list of all positions except the ones we want side chains for 
    # then iterate through list with mutate function

    # the pdb sequence for 3QN1 has gaps that we need to account for when using rosetta numbering 
    
    # niclosamide 
    # same_number_3QN1 = [60, 61, 62]
    # subtract2_3QN1 = [x - 2 for x in [79, 81, 83, 87, 88, 89, 90, 92, 94, 108, 109, 110, 115, 116, 117, 120, 122, 159, 160, 163, 164, 167]]
    # subtract26_3QN1 = [x - 26 for x in [385]]

    # WT_sidechain_positions = same_number_3QN1 + subtract2_3QN1 + subtract26_3QN1

    # mutant_sidechain_positions = [59, 139] # shifted 141 by 2 to account for gaps in sequence 
    # mutant_sidechain_dict = {59:'S', 139:'V'}

    # lynestrenol
    # 3QN1
    same_number_3QN1 = [59, 60, 61, 62]
    subtract2_3QN1 = [x - 2 for x in [79, 81, 83, 87, 88, 89, 90, 92, 94, 108, 109, 110, 115, 116, 117, 120, 122, 141, 160, 163, 167]]
    subtract26_3QN1 = [x - 26 for x in [385]]

    WT_sidechain_positions = same_number_3QN1 + subtract2_3QN1 + subtract26_3QN1

    # F159G V164S
    mutant_sidechain_positions = [157, 162] # shifted by 2 to account for gaps in sequence 
    mutant_sidechain_dict = {157:'G', 162:'S'}

    #4WVO
    # subtract2_4WVO = [x - 2 for x in [60, 61, 62]]
    # subtract6_4WVO = [x - 6 for x in [79, 81, 83, 87, 88, 89, 90, 92, 94, 108, 109, 110, 115, 116, 117, 120, 122, 141, 160, 163, 167]]
    # subtract31_4WVO = [x - 31 for x in [385]]

    # WT_sidechain_positions = subtract2_4WVO + subtract6_4WVO + subtract31_4WVO

    # # F159G V164S
    # mutant_sidechain_positions1 = [x - 6 for x in [159, 164]]
    # mutant_sidechain_positions2 = [x - 2 for x in [59]] 
    # mutant_sidechain_positions = mutant_sidechain_positions1 + mutant_sidechain_positions2
    # mutant_sidechain_dict = {57: 'G', 153:'G', 158:'S'}


    # alpha-estradiol
    # same_number_3QN1 = [59, 60, 61, 62]
    # subtract2_3QN1 = [x - 2 for x in [79, 81, 83, 87, 88, 89, 90, 92, 94, 108, 109, 110, 115, 116, 117, 120, 122, 141, 160, 167]]
    # subtract26_3QN1 = [x - 26 for x in [385]]

    # WT_sidechain_positions = same_number_3QN1 + subtract2_3QN1 + subtract26_3QN1

    # # F159V V163W V164G
    # mutant_sidechain_positions = [157, 161, 162] # shifted by 2 to account for gaps in sequence 
    # mutant_sidechain_dict = {157:'V', 161:'W', 162:'G'}

    # mandipropamid 4WVO
    # subtract2_4WVO = [x - 2 for x in [60, 61, 62]]
    # subtract6_4WVO = [x - 6 for x in [79, 83, 87, 88, 89, 90, 92, 94, 109, 110, 115, 116, 117, 120, 141, 160, 163, 164, 167]]
    # subtract31_4WVO = [x - 31 for x in [385]]

    # WT_sidechain_positions = subtract2_4WVO + subtract6_4WVO + subtract31_4WVO

    # # Y58H K59R V81I N90S F108A S122G F159L
    # mutant_sidechain_positions = [56, 57, 75, 84, 102, 116, 153]
    # mutant_sidechain_dict = {56:'H', 57:'R', 75:'I', 84:'S', 102:'A', 116:'G', 153:'L'}

    # isotonitazene
    # 3QN1
    same_number_3QN1 = [60, 61, 62]
    subtract2_3QN1 = [x - 2 for x in [79, 81, 83, 88, 89, 94, 109, 110, 115, 116, 117, 122, 141, 163, 164, 167]]
    subtract26_3QN1 = [x - 26 for x in [385]]

    WT_sidechain_positions = same_number_3QN1 + subtract2_3QN1 + subtract26_3QN1

    # #nitazene1 K59Q N90S S92M F108V Y120A F159A 
    # mutant_sidechain_positions = [59, 88, 90, 106, 118, 157] # shifted by 2 to account for gaps in sequence 
    # mutant_sidechain_dict = {59:'Q', 88:'S', 90:'M', 106:'V', 118:'A', 157:'A'}

    # #nitazene2 K59Q L87A N90S S92F F108V Y120A A160V V163F 
    # mutant_sidechain_positions = [59, 85, 88, 90, 106, 118, 158, 161] # shifted by 2 to account for gaps in sequence 
    # mutant_sidechain_dict = {59:'Q', 85:'A', 88:'S', 90:'F', 106:'V', 118:'A', 158:'V', 161:'F'}

    # #nitazene4 K59Q N90S S92L F108V Y120A V163F 
    # mutant_sidechain_positions = [59, 88, 90, 106, 118, 161] # shifted by 2 to account for gaps in sequence 
    # mutant_sidechain_dict = {59:'Q', 88:'S', 90:'L', 106:'V', 118:'A', 161:'F'}

    #nitazene8 K59Q L87A N90S S92M F108V Y120A F159A A160V 
    mutant_sidechain_positions = [59, 85, 88, 90, 106, 118, 157, 158] # shifted by 2 to account for gaps in sequence 
    mutant_sidechain_dict = {59:'Q', 85:'A', 88:'S', 90:'M', 106:'V', 118:'A', 157:'A', 158:'V'}


    # post_pose.dump_pdb('pre_mutate.pdb')

    # Iterate over each position in the peptide sequence
    for i in range(1, post_pose.total_residue() + 1):
        if i in WT_sidechain_positions:
            continue 
        if i in mutant_sidechain_positions:
            mutation = mutant_sidechain_dict[i]
            pyrosetta.toolbox.mutate_residue(post_pose, i, mutation)
        else: 
            pyrosetta.toolbox.mutate_residue(post_pose, i, 'G')

    post_pose.dump_pdb('post_mutate.pdb')


    # now back to regular code

    # Defining the residue to align to
    res = pre_pose.residue(pre_pose.pdb_info().pdb2pose(default["ChainLetter"], int(default["ResidueNumber"])))


    # Reading in information from the config file
    path_to_conformers = default["PathToConformers"]
    pkl_file_name = default["PKLFileName"]

    if pkl_file_name.strip() == "":
        pkl_file_name = None
    print(pkl_file_name)
    bin_width = float(spec["BinWidth"])
    vdw_modifier = float(spec["VDW_Modifier"])
    include_sc = spec["IncludeSC"] == "True"
    lig_res_num = int(default["LigandResidueNumber"])

    # Use an existant pkl file when possible
    print("Attempting to read in .pkl file")
    try:
        df = pd.read_pickle(pkl_file_name)
    except:
        print(".pkl file not found, generating one instead (this is normal)")
        csv_to_df_pkl(default["CSVFileName"], pkl_file_name, auto, path_to_conformers, pre_pose, res, lig_res_num)
        df = pd.read_pickle(pkl_file_name)

    print("\nBeginning grading")
    align_to_residue_and_check_collision(post_pose, res, path_to_conformers, df, pkl_file_name, 
                                        bin_width, vdw_modifier, include_sc, lig_res_num)

    

if __name__ == "__main__":  
    main(sys.argv[1:])
