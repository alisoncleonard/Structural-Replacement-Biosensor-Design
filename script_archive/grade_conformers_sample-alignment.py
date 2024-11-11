import sys
import os
import time

import pandas as pd
import numpy as np

from configparser import ConfigParser
import argparse
from math import radians

import pyrosetta
from pyrosetta.rosetta.core.pose import Pose
from pyrosetta.rosetta.core.select.residue_selector import ResidueIndexSelector
from pyrosetta.rosetta.core.id import AtomID
from pyrosetta.rosetta.numeric import xyzVector_double_t


# note for this application the rotation axis would need to be around the water - I think I have this for 3QN1 but can't remember how I calculated
# or the central target atom to the water? 

def rotate_ligand(pose, ligand_residue_ids, rotation_angle_degrees, axis):
    # Convert rotation angle to radians
    rotation_angle_rad = np.radians(rotation_angle_degrees)

    # Define rotation matrix
    axis = np.array(axis)
    axis = axis / np.linalg.norm(axis)
    c = np.cos(rotation_angle_rad)
    s = np.sin(rotation_angle_rad)
    t = 1 - c
    rotation_matrix = np.array([[t * axis[0]**2 + c, t * axis[0] * axis[1] - s * axis[2], t * axis[0] * axis[2] + s * axis[1]],
                                [t * axis[0] * axis[1] + s * axis[2], t * axis[1]**2 + c, t * axis[1] * axis[2] - s * axis[0]],
                                [t * axis[0] * axis[2] - s * axis[1], t * axis[1] * axis[2] + s * axis[0], t * axis[2]**2 + c]])

    # Apply rotation to the ligand atoms
    for residue_id in ligand_residue_ids:
        print(residue_id)
        for atom_id in range(1, pose.residue(residue_id).natoms() + 1):
            print(atom_id)
            atom_xyz = pose.residue(residue_id).xyz(atom_id)
            print(atom_xyz)
            atom_xyz_rotated = np.dot(rotation_matrix, atom_xyz)
            atom_xyz_rotated_doublet = xyzVector_double_t(atom_xyz_rotated[0], atom_xyz_rotated[1], atom_xyz_rotated[2])
            print(atom_xyz_rotated_doublet)
            pose.residue(residue_id).set_xyz(atom_id, atom_xyz_rotated_doublet)
            print(pose.residue(residue_id).xyz(atom_id))

def rotate_ligand_around_atom(pose, fixed_atom_id, ligand_residue_ids, rotation_angle_degrees, axis=[0, 0, 1]):
    # Get coordinates of the fixed atom
    fixed_atom_xyz = pose.residue(ligand_residue_ids[0]).xyz(fixed_atom_id)

    # Translate the ligand so that the fixed atom is at the origin
    for residue_id in ligand_residue_ids:
        for atom_id in range(1, pose.residue(residue_id).natoms() + 1):
            pose.residue(residue_id).set_xyz(atom_id, pose.residue(residue_id).xyz(atom_id) - fixed_atom_xyz)

    # Convert rotation angle to radians
    rotation_angle_rad = np.radians(rotation_angle_degrees)

    # Define rotation matrix
    axis = np.array(axis)
    axis = axis / np.linalg.norm(axis)
    c = np.cos(rotation_angle_rad)
    s = np.sin(rotation_angle_rad)
    t = 1 - c
    rotation_matrix = np.array([[t * axis[0]**2 + c, t * axis[0] * axis[1] - s * axis[2], t * axis[0] * axis[2] + s * axis[1]],
                                [t * axis[0] * axis[1] + s * axis[2], t * axis[1]**2 + c, t * axis[1] * axis[2] - s * axis[0]],
                                [t * axis[0] * axis[2] - s * axis[1], t * axis[1] * axis[2] + s * axis[0], t * axis[2]**2 + c]])

    # Apply rotation to the ligand
    for residue_id in ligand_residue_ids:
        for atom_id in range(1, pose.residue(residue_id).natoms() + 1):
            atom_xyz = pose.residue(residue_id).xyz(atom_id)
            atom_xyz_rotated = np.dot(rotation_matrix, atom_xyz)
            pose.residue(residue_id).set_xyz(atom_id, xyzVector_double_t(atom_xyz_rotated[0], atom_xyz_rotated[1], atom_xyz_rotated[2]))

    # Translate the ligand back to its original position
    for residue_id in ligand_residue_ids:
        for atom_id in range(1, pose.residue(residue_id).natoms() + 1):
            pose.residue(residue_id).set_xyz(atom_id, pose.residue(residue_id).xyz(atom_id) + fixed_atom_xyz)


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
         bin_width = 1, vdw_modifier = 0.7, include_sc = False, lig_res_num = 1, sample_angles = False):
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
    sample_angles: whether to sample variations on ligand alignment for each conformer

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

    axis_list = [[0, 1, 1], [0,0,1]]


    for pose_info in conformer_prep.yield_ligand_poses(df = df, path_to_conformers = path_to_conformers, post_accepted_conformers = False, ligand_residue = lig_res_num):
        
        if not pose_info:
            print(f"{conf.name}, {conf.id}: {len(accepted_conformations)/conf.conf_num}")
            all_accepted.append(accepted_conformations)
            accepted_conformations = []
            continue


        # Grab the conformer from the generator
        conf = pose_info

        
        # now that I have a way to rotate my conformer pose, I need to match the rotation with a loaded conformer
        # Need to make same number of rows for each conformer as I have angle arrays 

        # Perform alignment
        conf.align_to_target(res)

        # code to try and rotate ligand - need to update this eventually so that it only rotates when SampleAngles is false
        conf.pose.dump_pdb(f"aligned{total_confs}.pdb")

        rotation_angle_degrees = 30
        # Get residue IDs of ligand
        ligand_residue_ids = []
        for residue_id in range(1, conf.pose.total_residue() + 1):
            if not conf.pose.residue(residue_id).is_protein():
                ligand_residue_ids.append(residue_id)

        ligand_atom_types = []
        for residue_id in ligand_residue_ids: 
            residue = pose.residue(residue_id)
            # Iterate over atoms in the residue
            for atom_id in range(1, residue.natoms() + 1):
                atom_name = residue.atom_type(atom_id).element()
                # Add atom name to the list
                ligand_atom_types.append((residue_id, atom_id, atom_name))

        print(ligand_atom_types)


        # Rotate the ligand
        axis = axis_list[total_confs]
        # rotate_ligand(conf.pose, ligand_residue_ids, rotation_angle_degrees, axis)
        fixed_atom_id = 4
        rotate_ligand_around_atom(conf.pose, fixed_atom_id, ligand_residue_ids, rotation_angle_degrees, axis)
        conf.pose.dump_pdb(f"rotated{total_confs}.pdb")

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
    from pyrosetta.rosetta.core.select.residue_selector import ResidueIndexSelector
    from pyrosetta.rosetta.numeric import xyzVector_double_t
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
    #3QN1
    # same_number_3QN1 = [59, 60, 61, 62]
    # subtract2_3QN1 = [x - 2 for x in [79, 81, 83, 87, 88, 89, 90, 92, 94, 108, 109, 110, 115, 116, 117, 120, 122, 141, 160, 163, 167]]
    # subtract26_3QN1 = [x - 26 for x in [385]]

    # WT_sidechain_positions = same_number_3QN1 + subtract2_3QN1 + subtract26_3QN1

    # # F159G V164S
    # mutant_sidechain_positions = [157, 162] # shifted by 2 to account for gaps in sequence 
    # mutant_sidechain_dict = {157:'G', 162:'S'}

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


    # post_pose.dump_pdb('pre_mutate.pdb')

    # # Iterate over each position in the peptide sequence
    # for i in range(1, post_pose.total_residue() + 1):
    #     if i in WT_sidechain_positions:
    #         continue 
    #     if i in mutant_sidechain_positions:
    #         mutation = mutant_sidechain_dict[i]
    #         pyrosetta.toolbox.mutate_residue(post_pose, i, mutation)
    #     else: 
    #         pyrosetta.toolbox.mutate_residue(post_pose, i, 'G')

    # post_pose.dump_pdb('post_mutate.pdb')


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
    sample_angles = spec["SampleAngles"] == "False"

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
                                        bin_width, vdw_modifier, include_sc, lig_res_num, sample_angles)

    

if __name__ == "__main__":  
    main(sys.argv[1:])
