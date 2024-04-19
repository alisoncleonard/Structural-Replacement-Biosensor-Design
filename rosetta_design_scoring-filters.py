import sys
from configparser import ConfigParser
import argparse
import os

import pandas as pd
import numpy as np

import alignment
import conformer_prep
import collision_check


# Appends the ligand to the residue and sets constraints for it to not move as much
# Constraints kinda broke though

def prepare_pose_for_design(pose, focus_pos, lig, lig_aid):
    complx = pyrosetta.Pose()
    complx.assign(pose)

    vrm = pyrosetta.rosetta.protocols.simple_moves.VirtualRootMover()
    vrm.apply(complx)

    # complx.append_residue_by_jump(lig.residue(1), focus_pos, "", "", True)
    complx.append_residue_by_jump(lig.residue(1), focus_pos, "", "", False)
    # complx.append_residue_by_jump(lig.residue(1), focus_pos, "CE2", "O2", False)

    lig_pos = pyrosetta.rosetta.core.select.residue_selector.ResidueIndexSelector()
    lig_pos.set_index(complx.size())

    complx.pdb_info().chain(complx.size(), 'X')
    complx.pdb_info().number(complx.size(), 1)
    complx.pdb_info().obsolete(False)

    ccg = pyrosetta.rosetta.protocols.constraint_generator.CoordinateConstraintGenerator()
    ccg.set_residue_selector(lig_pos)
    ccg.set_sidechain(True)
    ccg.set_bounded_width(0.2)
    ccg.set_sd(0.25)
    ccg.apply(complx)

    return complx

# Determine the mutations from a string of mutations and return them
# IE something like W100G A502V
def determine_mutations(pose, mutation_string):
    if mutation_string == "WT":
        pos_res = []
    else:
        mutations = mutation_string.split(" ")
        pos_res = [(int(mutation[1:-1]), mutation[-1:]) for mutation in mutations]

    return [(pose.pdb_info().pdb2pose("A", pos), res) for pos, res in pos_res]

def design(pose, pose_nolig, mt_pose_nolig, wt_pose_ligand, focus_seqpos, resfile, ref_seq, conf_id, scoring_filters_file, fnr_bonus):

    # SET UP SCORE FUNCTION
    sf = pyrosetta.get_fa_scorefxn()
    sf.set_weight(pyrosetta.rosetta.core.scoring.res_type_constraint, 1.0)
    sf.set_weight(pyrosetta.rosetta.core.scoring.coordinate_constraint, 1.0)
    print(sf)

    # score without buried unsatisfied penalty to make cluster work
    # unsat_penalty = 1.0
    # sf.set_weight(pyrosetta.rosetta.core.scoring.buried_unsatisfied_penalty, unsat_penalty)

    # SET UP SCORE BONUSES
    fnr_bonus = 4
    fnr = pyrosetta.rosetta.protocols.protein_interface_design.FavorNativeResidue(pose, fnr_bonus)
    print('FNR bonus: ' + str(fnr_bonus))


    # seq_profile = pyrosetta.rosetta.core.sequence.SequenceProfile()
    # pssm_file = '../../pssm/test_PYR1-identity-short-logodds.pssm'
    # print(pssm_file)
    # seq_profile.read_from_file(pyrosetta.rosetta.utility.file.FileName(pssm_file))
    # fsp = pyrosetta.rosetta.protocols.simple_moves.FavorSequenceProfile()
    # fsp.set_profile(seq_profile)
    # fsp_weight = 1.0
    # fsp.set_weight(fsp_weight)
    # print('fsp weight: ' + str(fsp_weight))
    # fsp.apply(pose)


    # SET UP TASK FACTORY
    tf = pyrosetta.rosetta.core.pack.task.TaskFactory()
    tf.push_back(pyrosetta.rosetta.core.pack.task.operation.InitializeFromCommandline())
    tf.push_back(pyrosetta.rosetta.core.pack.task.operation.IncludeCurrent())
    tf.push_back(pyrosetta.rosetta.core.pack.task.operation.NoRepackDisulfides())
    tf.push_back(pyrosetta.rosetta.core.pack.task.operation.ReadResfile(resfile))

    # Do not use neighborhood residues, code is restricting design to only these residues incorrectly
    # SET UP NEIGHBORHOOD RESIDUES and JUMP
    # focus_pos = pyrosetta.rosetta.core.select.residue_selector.ResidueIndexSelector()
    # focus_pos.set_index(focus_seqpos)
    # nbr = pyrosetta.rosetta.core.select.residue_selector.NeighborhoodResidueSelector()
    # nbr.set_distance(12)
    # # nbr.set_distance(100)
    # nbr.set_focus_selector(focus_pos)
    # nbr.set_include_focus_in_subset(True)
    #
    # pose.update_residue_neighbors()
    # nbr.apply(pose)
    # rlt = pyrosetta.rosetta.core.pack.task.operation.PreventRepackingRLT() # RLT = residue level task operators
    # tf.push_back(pyrosetta.rosetta.core.pack.task.operation.OperateOnResidueSubset(rlt, nbr, True))

    jump_num = pose.num_jump()
    jump_selector = pyrosetta.rosetta.core.select.jump_selector.JumpIndexSelector()
    jump_selector.jump(jump_num)
    print('number jumps', jump_num)


    # SET UP MOVEMAP
    mmf = pyrosetta.rosetta.core.select.movemap.MoveMapFactory()
    # mmf.add_chi_action(pyrosetta.rosetta.core.select.movemap.mm_enable, nbr)
    mmf.add_jump_action(pyrosetta.rosetta.core.select.movemap.mm_enable, jump_selector)

    print('pre-design score: ' + str(sf(pose)))

    # SET UP FASTRELAX
    fr = pyrosetta.rosetta.protocols.relax.FastRelax(sf, standard_repeats = 3)
    fr.set_task_factory(tf)
    fr.set_movemap_factory(mmf)
    fr.ramp_down_constraints(False)
    fr.min_type('lbfgs_armijo_nonmonotone')
    fr.max_iter(300)

    print('finish setup')

    print(pose.sequence())
    print(pose.fold_tree())

    #Apply FastRelax
    fr.apply(pose)

    print('finish apply fastrelax')
    print('post-design score: ' + str(sf(pose)))

    # Minimization without the buried unsatisfied penalty, recc. by Vmulligan

    sf.set_weight(pyrosetta.rosetta.core.scoring.res_type_constraint, 0)

    repack_tf = pyrosetta.rosetta.core.pack.task.TaskFactory()
    repack_tf.push_back(pyrosetta.rosetta.core.pack.task.operation.InitializeFromCommandline())
    repack_tf.push_back(pyrosetta.rosetta.core.pack.task.operation.IncludeCurrent())
    repack_tf.push_back(pyrosetta.rosetta.core.pack.task.operation.NoRepackDisulfides())
    repack_tf.push_back(pyrosetta.rosetta.core.pack.task.operation.RestrictToRepacking())
    # repack_tf.push_back(pyrosetta.rosetta.core.pack.task.operation.OperateOnResidueSubset(rlt, nbr, True))

    fr = pyrosetta.rosetta.protocols.relax.FastRelax(sf, standard_repeats = 3)
    fr.set_task_factory(repack_tf)
    fr.set_movemap_factory(mmf)
    fr.ramp_down_constraints(False)
    fr.min_type('lbfgs_armijo_nonmonotone')
    fr.max_iter(300)

    # sf.set_weight(pyrosetta.rosetta.core.scoring.buried_unsatisfied_penalty, 0)

    fr.set_scorefxn(sf)
    fr.apply(pose)

    print('finish Minimization')
    post_min_score = sf(pose)
    print('post-minimization score: ' + str(post_min_score))



    # Determine mutations
    # this code is just giving you the printout of where the suggested mutations are, and how they different from vdm mutations
    new_seq = pose.sequence()
    pdb_num_mutations = []
    rosetta_mutations = []
    for i in range(len(ref_seq)): # ref_seq is the AA sequence from DesignPDBFile
        if new_seq[i] != ref_seq[i]:
            empty_string = " "
            pdb_num_mutations.append(f"{ref_seq[i]}{pose.pdb_info().pose2pdb(i+1).split(empty_string)[0]}({i+1}){new_seq[i]}")
            rosetta_mutations.append(f"{ref_seq[i]}{pose.pdb_info().pose2pdb(i+1).split(empty_string)[0]}{new_seq[i]}")

    rosetta_mutations = " ".join(rosetta_mutations)


    print('start interface scoring filters')
    # https://www.rosettacommons.org/docs/latest/scripting_documentation/RosettaScripts/Filters/Filters-RosettaScripts#special-application-filters
    # https://www.rosettacommons.org/docs/latest/scripting_documentation/RosettaScripts/Movers/movers_pages/InterfaceScoreCalculatorMover
    # https://graylab.jhu.edu/PyRosetta.documentation/pyrosetta.rosetta.protocols.scoring.html#pyrosetta.rosetta.protocols.scoring.Interface
    # https://github.com/RosettaCommons/PyRosetta.notebooks/blob/master/notebooks/08.01-Ligand-Docking-XMLObjects.ipynb

    xml = pyrosetta.rosetta.protocols.rosetta_scripts.XmlObjects.create_from_file(scoring_filters_file).get_mover("ParsedProtocol")
    xml.apply(pose)
    print(pose.scores)

    interfE = pose.scores['interfE']
    sasa = pose.scores['sasa']
    buried_unsats_ligjump = pose.scores['buried_unsat_Hbonds_ligjump']
    ligand_hbonds = pose.scores['hbonds_to_ligand']
    shape_comp = pose.scores['ShapeComp']
    total_score = pose.scores['total_score']


    return (pose, rosetta_mutations, pdb_num_mutations, interfE, sasa, buried_unsats_ligjump, ligand_hbonds, shape_comp, total_score)



def main(argv):
    # PARSING CONFIG FILE
    print("Parsing config file")

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("config_file",
                        help = "your config file",
                        default = "my_conf.txt")

    if len(argv) == 0:
        print(parser.print_help())
        return

    args = parser.parse_args(argv)

    config = ConfigParser()
    config.read(args.config_file)
    default = config["DEFAULT"]
    spec = config["rosetta_design"]


    # Importing necessary dependencies
    print("Importing necessary dependencies")
    sys.path.append(default["PathToPyRosetta"])
    global pyrosetta, Pose
    import pyrosetta
    from pyrosetta.rosetta.core.pose import Pose

    pyrosetta.init("-mute all -multithreading:total_threads 10")

    # Grabbing values from config file
    print("Grabbing info from config file")
    params_list = default["ParamsList"]
    path_to_conformers = default["PathToConformers"]
    path_to_complexes = spec["PathToComplexes"]
    scoring_filters_file = spec["Filters"]
    subset_start = int(spec['SubStart'])
    subset_end = int(spec['SubEnd'])


    # make design output folder
    try:
        print(f"Attempting to make {path_to_complexes}")
        os.mkdir(path_to_complexes)
    except:
        print(f"{path_to_complexes} already made, continuing")

    resfile = spec["Resfile"]
    fnr_bonus = float(spec["FNRBonus"])
    vdm_bonus = float(spec["VDMBonus"])
    # the vdm output file you wish to design against
    pkl_file = spec["DesignPKLFile"]

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

    # structural feature to map on to
    res = pre_pose.residue(pre_pose.pdb_info().pdb2pose(default["ChainLetter"], int(default["ResidueNumber"])))

    # "focus" residue identified in config file
    focus_pos = pre_pose.pdb_info().pdb2pose(spec["FocusChain"], int(spec["FocusResidueNumber"]))

    new_df = []
    df = pd.read_pickle(pkl_file)

    print("\n\nBeginning Design\n\n")
    score_type = "Norm. Score"
    df = df.sort_values(score_type, ascending = False)#.drop_duplicates(["Molecule ID"], keep = "first")
    # note: no sorting if there are no scores
    # removed drop_duplicates['Molecule ID'] because this doesn't work when I'm testing lots of different alignments of the same conformers

    # for i, molecule_df in df.iterrows():
    for i, molecule_df in df[subset_start:subset_end].iterrows():
        # can change range of row depending on what you want to design around. Future task will be selecting which conformers from list to design, particularly when there is no vdm score
        # note: design only happens on the first # of specified conformers, so if there are more conformers but the conformers have no score (because no mutations) then you will score a potentially unrepresentative subset

        print("\n\nPreparing pose for design")
        pose = pyrosetta.pose_from_pdb(spec["DesignPDBFile"])
        pose_nolig = pyrosetta.pose_from_pdb(spec["DesignPDBFile"])
        mt_pose_nolig = pyrosetta.pose_from_pdb(spec["DesignPDBFile"])
        ref_seq = pose.sequence()

        mol_id, conf_num, file_stem, lig_aid, r_aid, mutation_string, vdm_positions = molecule_df[[
        "Molecule ID", "Conformer Number", "Molecule File Stem", "Molecule Atoms", "Target Atoms",
        "Mutations", "Positions"]]

        # pose_nolig.dump_pdb(f"{path_to_complexes}/dump_noligand.pdb")

        print(f"{mol_id}, {conf_num}")

        lig = Pose()
        res_set = pyrosetta.generate_nonstandard_residue_set(lig, params_list = [f"{path_to_conformers}/{file_stem}.params"])
        pyrosetta.pose_from_file(lig, res_set, f"{path_to_conformers}/{file_stem}_{conf_num:04}.pdb")

        alignment.align_ligand_to_residue(lig, res, lig_aid, r_aid)

        pose_with_ligand = prepare_pose_for_design(pose, focus_pos, lig, lig_aid)
        wt_pose_ligand = prepare_pose_for_design(pose, focus_pos, lig, lig_aid)

        conf_id = f"{mol_id}_{conf_num}"
        #pose_with_ligand.dump_pdb(f"{path_to_complexes}/{mol_id}_{conf_num:04}_angle{angle}_translate{t}_dump_predesign.pdb")

        # update this to include WT pose w/o ligand, need for scoring ddg. Then can make mutations to w/o ligand scaffold using "mutate_residue"
        out_pose, rosetta_mutations, pdb_num_mutations, interfE, sasa, buried_unsats_ligjump, ligand_hbonds, complementarity, total_score = design(pose_with_ligand, pose_nolig, mt_pose_nolig, wt_pose_ligand, focus_pos, resfile, ref_seq, conf_id, scoring_filters_file, fnr_bonus = fnr_bonus)

        new_df.append(list(molecule_df[:-2]) + [pdb_num_mutations, rosetta_mutations, interfE, sasa, buried_unsats_ligjump, ligand_hbonds, complementarity, total_score])
        # out_pose.dump_pdb(f"{path_to_complexes}/{mol_id}_{conf_num:04}_angle{angle}.pdb")
        out_pose.dump_pdb(f"{path_to_complexes}/{mol_id}_{conf_num:04}.pdb")


    new_df = pd.DataFrame(new_df)

    new_df.columns = list(df.columns)[:-2] + ["Pose Numbered Mutations", "Rosetta Mutations", "interface energy", 'delta SASA', "buried unsat Hbonds ligjump", "Hbonds to ligand", "complementarity", "total score"]
    new_df.to_pickle(f"{path_to_complexes}/complex_scores.pkl")




if __name__ == "__main__":
    main(sys.argv[1:])
