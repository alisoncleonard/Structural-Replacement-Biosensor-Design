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

def align_to_residue_and_check_collision(pose, res, path_to_conformers, df, pkl_file, jump_num, rotation, translation, upper_water_distance, 
         lower_water_distance, backbone_clash_keep_sidechains, max_pass_score = -300, bin_width = 1, vdw_modifier = 0.7, include_sc = False, lig_res_num = 1):
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

    import pyrosetta.rosetta.protocols.rigid as rigid_moves
    import pyrosetta.rosetta.protocols.grafting as graft
    from pyrosetta.rosetta.core.pack.task import TaskFactory
    from pyrosetta.rosetta.core.pack.task import operation
    from pyrosetta.rosetta.protocols import minimization_packing as pack_min

    pmm = pyrosetta.PyMOLMover()
    pmm.keep_history(True)

    # Setup all scorefunction and packing stuff here outside loop
    tf = TaskFactory()

    tf.push_back(operation.InitializeFromCommandline())
    tf.push_back(operation.IncludeCurrent())
    tf.push_back(operation.NoRepackDisulfides())
    tf.push_back(operation.RestrictToRepacking())
    packer = pack_min.PackRotamersMover()
    packer.task_factory(tf)
    
    sf_all = pyrosetta.get_fa_scorefxn()
    sf_rep = pyrosetta.get_fa_scorefxn()
    sf_rep.set_weight(pyrosetta.rosetta.core.scoring.ScoreType.fa_rep, 1.0)  # Set weight for the van der Waals repulsive term

    t0 = time.time()

    all_accepted = []
    accepted_conformations = []
    total_confs = 0

    # Set up post_pose for merged PDB file
    alignto_pose = pose

    # Make a new pose for backbone grid check that includes side chains on some loops.
    loop_include_sidechain_pose = pose.clone()
    loop_include_sidechain = [int(x) for x in backbone_clash_keep_sidechains]

    # Iterate over each position in the peptide sequence
    for i in range(1, loop_include_sidechain_pose.total_residue() + 1):
        if i in loop_include_sidechain:
            continue 
        else: 
            pyrosetta.toolbox.mutate_residue(loop_include_sidechain_pose, i, 'G')

    # Precalculating the protein collision grid and storing the results in a hashmap (grid.grid)
    backbone_grid = collision_check.CollisionGrid(loop_include_sidechain_pose, bin_width = bin_width, vdw_modifier = vdw_modifier, include_sc = True)

    # Iterate over all conformers specified in Ligands.csv, must make as many conformers as I want to iterate through
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

        # Create merged PDB file for ligand perturbation
        location = len(alignto_pose.chain_sequence(1))  # getting lenth of first chain in pose 
        new_pose = graft.insert_pose_into_pose(alignto_pose, conf.pose, location)

        # print(new_pose.fold_tree())

        ligand_res_index = location + 1

        # Try making ligand perturbations until I have one that satisfies criteria 

        keep_list = []
        count = 0
        while len(keep_list) == 0 and count < 10:  # pick first conformer perturbation that satisfies criteria but stop after 10 tries if nothing works 


            copy_pose = new_pose.clone() # need to copy pose each time or else the next perturbations use the previous as the starting point 

            # Do ligand perturbation in pocket 
            pert_mover = rigid_moves.RigidBodyPerturbMover(jump_num, rotation, translation)
            pert_mover.apply(copy_pose)
            
            # Check if the ligand collides with the backbone only 

            # Extract the new ligand coordinates and update coordinates of conf
            for atom_id in range(1, copy_pose.residue(ligand_res_index).natoms() + 1):
                atom_coords = copy_pose.residue(ligand_res_index).xyz(atom_id)
                conf.pose.residue(1).set_xyz(atom_id, atom_coords)

            does_collide = conf.check_collision(backbone_grid) 

            if not does_collide:

                # Check if new ligand alignment is close enough to water to make H-bond

                water_resid = []
                # need to get coordinates of TP3 residues, keep them in a list that I could iterate through to accomodate single or multiple waters
                for i in range(1, copy_pose.total_residue()+1):
                    residue = copy_pose.residue(i)
                    if residue.is_water():
                        water_resid.append(i)

                # include if statement to skip this step if an H-bond to water is not required for binding 
                if len(water_resid) == 0:
                    pass
                else: 
                    for j in range(1, copy_pose.residue(ligand_res_index).natoms() + 1):
                        atom = copy_pose.residue(ligand_res_index).atom_type(j)

                        if atom.element() == "N" or atom.element() == "O":
                            # print(atom.element())
                            atom_coords = copy_pose.residue(ligand_res_index).xyz(j)
                    
                            for i in range(len(water_resid)): 
                                # print(water_resid[i])
                                for w in range(1, copy_pose.residue(water_resid[i]).natoms() + 1):
                                    atom = copy_pose.residue(water_resid[i]).atom_type(w)
                                    if atom.element() == "O":
                                        water_coord = copy_pose.residue(water_resid[i]).xyz(w)
                                distance_to_water = atom_coords.distance(water_coord)
                                # print(distance_to_water)

                                # 2.7 to 3.3 angstroms seemed to miss some alignments in PyMol 
                                if distance_to_water < upper_water_distance and distance_to_water > lower_water_distance:
                                    keep_list.append(i)
            # print(count)
            count += 1

        if len(keep_list) == 0: 
            print("conformer failed backbone clash and water distance check")
        else: 

            copy_pose.dump_pdb(f"rotated/rotated{total_confs}.pdb")
            
            # Repack and minimize 

            packer.apply(copy_pose)

            # copy_pose.dump_pdb(f"repacked/repacked{total_confs}.pdb")

            # Extract the new ligand coordinates and update coordinates of conf
            ligand_res_index = location + 1
            for atom_id in range(1, copy_pose.residue(ligand_res_index).natoms() + 1):
                atom_coords = copy_pose.residue(ligand_res_index).xyz(atom_id)
                conf.pose.residue(1).set_xyz(atom_id, atom_coords)

            # conf.pose.dump_pdb(f"conf{total_confs}.pdb")

            # Check for collision
            remove_ligand_pose = copy_pose.clone()
            remove_ligand_pose.delete_residue_slow(ligand_res_index)
            grid = collision_check.CollisionGrid(remove_ligand_pose, bin_width = bin_width, vdw_modifier = vdw_modifier, include_sc = include_sc)
            does_collide = conf.check_collision(grid) 

            if not does_collide:
                accepted_conformations.append(conf.conf_num)

                score = sf_all(copy_pose)
                print(f"molecule repacked {total_confs}", score)

                if score < max_pass_score:

                    # View alignment for conformers passing collision check
                    conf.pose.pdb_info().name(f"pass_score")
                    pmm.apply(conf.pose)
                    copy_pose.dump_pdb(f"pass_score_repacked/repacked{total_confs}.pdb")

                else: 
                    # View alignment for conformers passing collision check
                    conf.pose.pdb_info().name(f"fail_score")
                    pmm.apply(conf.pose)
                    copy_pose.dump_pdb(f"pass_fail-score_repacked/repacked{total_confs}.pdb")
            else: 
                copy_pose.dump_pdb(f"pass_initial-only_repacked/repacked{total_confs}.pdb")

        total_confs += 1

        

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
    parser.add_argument("-c", "--passes_collision_check",
                        help = "Whether to only show ones that pass collison check (True/False)",
                        default = "False")
    
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
    global pyrosetta, Pose, alignment, conformer_prep, collision_check, csv_to_df_pkl

    import pyrosetta
    from pyrosetta.rosetta.core.pose import Pose
    import alignment
    import conformer_prep
    import collision_check
    from grade_conformers_docked_to_sequence import csv_to_df_pkl

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

    
    # Set protein sequence for docking  

    positions = str(spec["MutatePositions"]).split()
    # convert a list of strings into a list of integers 
    mutant_sidechain_positions = [int(x) for x in positions]

    AA_values = str(spec["MutateAA"]).split()
    mutant_sidechain_dict = dict(zip(mutant_sidechain_positions, AA_values))

    # just mutations no glycine shave
    for i in range(1, post_pose.total_residue() + 1):
        if i in mutant_sidechain_positions:
            mutation = mutant_sidechain_dict[i]
            pyrosetta.toolbox.mutate_residue(post_pose, i, mutation)

    post_pose.dump_pdb('post_mutate.pdb')

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
    # Info for ligand perturbation 
    jump_num = int(spec["JumpNum"])
    rotation = float(spec["Rotation"])        
    translation = float(spec["Translation"])  
    upper_water_distance = float(spec["UpperWaterDistance"])
    lower_water_distance = float(spec["LowerWaterDistance"])
    backbone_clash_keep_sidechains = str(spec["BackboneClashKeepSidechains"]).split()
    max_pass_score = float(spec["MaxScore"])

    # Use an existant pkl file when possible
    print("Attempting to read in .pkl file")
    try:
        df = pd.read_pickle(pkl_file_name)
    except:
        print(".pkl file not found, generating one instead (this is normal)")
        csv_to_df_pkl(default["CSVFileName"], pkl_file_name, auto, path_to_conformers, pre_pose, res, lig_res_num)
        df = pd.read_pickle(pkl_file_name)

    try:
        os.mkdir('./rotated')
    except:
        print("didn't make directory rotated")

    try:
        os.mkdir('./pass_score_repacked')
    except:
        print("didn't make directory pass_score_repacked")

    try:
        os.mkdir('./pass_fail-score_repacked')
    except:
        print("didn't make directory pass_fail-score_repacked")

    try:
        os.mkdir('./pass_initial-only_repacked')
    except:
        print("didn't make directory pass_initial-only_repacked")


    print("\nBeginning grading")
    align_to_residue_and_check_collision(post_pose, res, path_to_conformers, df, pkl_file_name, jump_num, rotation, translation, upper_water_distance, 
                                        lower_water_distance, backbone_clash_keep_sidechains, max_pass_score, bin_width, vdw_modifier, include_sc, lig_res_num)

    

if __name__ == "__main__":  
    main(sys.argv[1:])
