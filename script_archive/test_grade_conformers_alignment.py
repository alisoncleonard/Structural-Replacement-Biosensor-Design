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

    alignto_pose = pose

    # print(alignto_pose.fold_tree())
    # total_residues = alignto_pose.total_residues()
    # print(total_residues)


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

        # code to try and rotate ligand - need to update this eventually so that it only rotates when SampleAngles is false
        conf.pose.dump_pdb(f"aligned{total_confs}.pdb")

        alignto_pose.dump_pdb('pyr1_pose.pdb')

        new_pose = pyrosetta.rosetta.protocols.grafting.insert_pose_into_pose(alignto_pose, conf.pose, 475)

        print(new_pose.fold_tree())

        # I should make sure that this code works for rotation and translation before I try and implement it here



        # rotation_angle_degrees = 30
        # # Get residue IDs of ligand
        # ligand_residue_ids = []
        # for residue_id in range(1, conf.pose.total_residue() + 1):
        #     if not conf.pose.residue(residue_id).is_protein():
        #         ligand_residue_ids.append(residue_id)

        # ligand_atom_types = []
        # for residue_id in ligand_residue_ids: 
        #     residue = pose.residue(residue_id)
        #     # Iterate over atoms in the residue
        #     for atom_id in range(1, residue.natoms() + 1):
        #         atom_name = residue.atom_type(atom_id).element()
        #         # Add atom name to the list
        #         ligand_atom_types.append((residue_id, atom_id, atom_name))

        # print(ligand_atom_types)


        # # Rotate the ligand
        # axis = axis_list[total_confs]
        # # rotate_ligand(conf.pose, ligand_residue_ids, rotation_angle_degrees, axis)
        # fixed_atom_id = 4
        # rotate_ligand_around_atom(conf.pose, fixed_atom_id, ligand_residue_ids, rotation_angle_degrees, axis)
        # conf.pose.dump_pdb(f"rotated{total_confs}.pdb")

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

def test_alignment(pose, res, path_to_conformers, df, pkl_file,
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
    import pyrosetta.rosetta.protocols.rigid as rigid_moves 


    t0 = time.time()

    # used old docking tutorial to implement RigidBodyPerturbMover 
    # https://graylab.jhu.edu/pyrosetta/downloads/documentation/Workshop7_PyRosetta_Docking.pdf

    ligand_pose = pose

    print(ligand_pose.fold_tree())
    ligand_pose.dump_pdb('pyr1_pose.pdb')
    # residue 181 is the water, so 182 must be the ligand 

    for count in range(10):


    # pyrosetta.rosetta.protocols.docking.setup_foldtree(ligand_pose, "ABD_C", Vector1([1]))
    # print(ligand_pose.fold_tree())

    # print(ligand_pose.jump(3).get_rotation())
    # print(ligand_pose.jump(3).get_translation())

        jump_num = 1
        rotation = 16     # mean rotation angle in degrees 
        translation = 1.2  # mean translation in angstroms
        pert_mover = rigid_moves.RigidBodyPerturbMover(jump_num, rotation, translation)

        ligand = ligand_pose.clone()
        pert_mover.apply(ligand)

        ligand.dump_pdb(f"dump{count}.pdb")



    # for pose_info in conformer_prep.yield_ligand_poses(df = df, path_to_conformers = path_to_conformers, post_accepted_conformers = False, ligand_residue = lig_res_num):
        
    #     if not pose_info:
    #         print(f"{conf.name}, {conf.id}: {len(accepted_conformations)/conf.conf_num}")
    #         all_accepted.append(accepted_conformations)
    #         accepted_conformations = []
    #         continue

    #     # Grab the conformer from the generator
    #     conf = pose_info
        
    #     # Perform alignment
    #     conf.align_to_target(res)

    #     # code to try and rotate ligand - need to update this eventually so that it only rotates when SampleAngles is false
    #     conf.pose.dump_pdb(f"aligned{total_confs}.pdb")

    #     alignto_pose.dump_pdb('pyr1_pose.pdb')

    #     new_pose = pyrosetta.rosetta.protocols.grafting.insert_pose_into_pose(alignto_pose, conf.pose, 475)

    #     print(new_pose.fold_tree())

        # I should make sure that this code works for rotation and translation before I try and implement it here




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
    # align_to_residue_and_check_collision(post_pose, res, path_to_conformers, df, pkl_file_name, 
                                        # bin_width, vdw_modifier, include_sc, lig_res_num, sample_angles)
    test_alignment(pre_pose, res, path_to_conformers, df, pkl_file_name, 
                                        bin_width, vdw_modifier, include_sc, lig_res_num, sample_angles)

    

if __name__ == "__main__":  
    main(sys.argv[1:])
