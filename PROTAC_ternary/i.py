import numpy as np

from Args import args


import scipy.optimize
import math


import itertools






def atom_list_generation(atom_file):
    """Take in atom files and output atom_names and residue_IDs list, for both 
       alignment atom files and delete atom files.""" 
    atom_list = [line.rstrip('\n').split() for line in open(atom_file)]
    return atom_list


def pdb_to_dictionary(pdb_file):
    """Take in linker/decoy pdb and output a dictionary with (ligand name + 
       atom name) as key and (x, y, z) as val."""

    pdb_dic = {}

    with open(pdb_file, "r") as f_open:
        for line in f_open:
            if line.startswith("HETATM"):
                res_name = line[17:20].strip()
                atom_name = line[12:16].strip()
                x_coor = float(line[30:38].strip())
                y_coor = float(line[38:46].strip())
                z_coor = float(line[46:54].strip())
                pdb_dic[res_name + " " + atom_name] = [x_coor, y_coor, z_coor]

    return pdb_dic


def pdb_coors_assignment(pdb_list, pdb_dic):
    """Take in dictionary and output coordinates array (6X3)."""
    pdb_coors = np.empty([6, 3])


    for i, atom in enumerate(pdb_list):
        x_coor, y_coor, z_coor = [],[],[]

        ligand_name = atom[0]

        for atom_name in atom[1::]:
            coors = pdb_dic[ligand_name + " " + atom_name]
            x_coor.append(coors[0])
            y_coor.append(coors[1])
            z_coor.append(coors[2])

        pdb_coors[i, 0] = np.average(x_coor)
        pdb_coors[i, 1] = np.average(y_coor)
        pdb_coors[i, 2] = np.average(z_coor)

    return pdb_coors












def apply_transformation(xtrans_ytrans_ztrans_rot1_rot2_rot3, orig_coors):
    """Take in original coordinates and transformation vetors, and output new coordinates after transformation."""
    # transformation will consist of first rotating, then translating
    new_coors = np.copy(orig_coors)
    xtrans, ytrans, ztrans, rot1, rot2, rot3 = xtrans_ytrans_ztrans_rot1_rot2_rot3
    # consolidate the three rotations into a single transformation matrix
    rot_mat1 = np.array([[1, 0, 0], [0, math.cos(
        rot1), -math.sin(rot1)], [0, math.sin(rot1), math.cos(rot1)]])
    rot_mat2 = np.array([[math.cos(rot2), 0, math.sin(rot2)], [
                        0, 1, 0], [-math.sin(rot2), 0, math.cos(rot2)]])
    rot_mat3 = np.array([[math.cos(rot3), -math.sin(rot3), 0],
                         [math.sin(rot3), math.cos(rot3), 0], [0, 0, 1]])
    rot_mat_all = np.dot(np.dot(rot_mat1, rot_mat2), rot_mat3)
    # apply the rotation matrix
    for index in range(len(new_coors[:, 0])):
        new_coors[index, :] = rot_mat_all.dot(new_coors[index, :])
    # apply the translation
    new_coors[:, 0] += xtrans
    new_coors[:, 1] += ytrans
    new_coors[:, 2] += ztrans
    return new_coors

def eval_rmsd_after_transformation(xtrans_ytrans_ztrans_rot1_rot2_rot3, moving_coors, ref_coors):
    """Take in two groups of coordinates, do transformation for one of them, and calculate the sum of the squares of the error as output."""
    xtrans, ytrans, ztrans, rot1, rot2, rot3 = xtrans_ytrans_ztrans_rot1_rot2_rot3
    new_coors = apply_transformation(
        xtrans_ytrans_ztrans_rot1_rot2_rot3, moving_coors)
    dist = new_coors - ref_coors
    dist = dist * dist
    sq_err = np.sum(dist)
    # note: this is not actually the RMSD, it's the sum of the squares of the error (for speed)
    return sq_err


if __name__ == "__main__":
    args = args()

    if args.cutoff < 0:
        print("Error: cut_off must be positive.")
        quit()


    decoys_list  = args.decoy
    linkers_list = args.linker
    cutoff       = args.cutoff
    rmsd_file    = args.rmsd

    linker_alignment_list = atom_list_generation(args.linker_alignment)
    decoy_alignment_list  = atom_list_generation(args.decoy_alignment)

    warheads_atom_delete_list  = atom_list_generation(args.warheads_delete)
    linker_atom_delete_list    = atom_list_generation(args.linker_delete)

    
    decoy_coords = {}

    for fname in decoys_list:
        dic = pdb_to_dictionary(fname)
        coords = pdb_coors_assignment(decoy_alignment_list, dic)
        decoy_coords[fname] = coords
    # print(str(len(decoy_list)) + " decoy(s) has/have been passed to the system.")

    linker_coords = {}

    for fname in linkers_list:
        dic = pdb_to_dictionary(fname)
        coords = pdb_coors_assignment(linker_alignment_list, dic)
        linker_coords[fname] = coords
    # print(str(len(linker_list)) + " linker conformer(s) has/have been passed to the system.")

    list_to_write = []
    ternary_model_number = 0

    for (d_name, l_name) in itertools.product(decoy_coords, linker_coords):
        d_coords = decoy_coords[d_name]
        l_coords = linker_coords[l_name]
        natoms = len(l_coords[:, 0])


        res = scipy.optimize.minimize(eval_rmsd_after_transformation,
                                     (0., 0., 0., 0., 0., 0.),
                                     args=(l_coords, d_coords),
                                     method='Powell')

        new_l_coords = apply_transformation((res.x[0], res.x[1], res.x[2], res.x[3], res.x[4], res.x[5]), l_coords)

        
        rmsd_value = eval_rmsd_after_transformation((0., 0., 0., 0., 0., 0.),
                                                    new_l_coords,
                                                    d_coords)
        
        rmsd_value = math.sqrt(rmsd_value/natoms)
        print(rmsd_value)
        break


    #     write_line = f"{d_name} {l_name} {str(rmsd_value)}"

    #     list_to_write.append(write_line + "\n")



    #     # print(rmsd_value, cutoff)
    #     # if rmsd_value < cutoff and ternary_model_name != "":



    #     #     if ternary_model_name == "default":
    #     #         ternary_model = "ternary" + str(ternary_model_number) + ".pdb"
    #     #     else:
    #     #         ternary_model = (decoy_pdb.split("/")[-1]).split(".")[0] + "_" + linker_conformer_pdb.split("/")[-1]

    #     #     line_number = len(linker_wholefile_dic)
    #     #     orig_whole_linker_coors = whole_linker_coors_assignment(linker_wholefile_dic, line_number)
    #     #     new_whole_linker_coors = apply_transformation((res.x[0], res.x[1], res.x[2], res.x[3], res.x[4], res.x[5]), orig_whole_linker_coors)
    #     #     new_linker_wholefile_dic = linker_conformer_after_rmsd_optimizaton(linker_wholefile_dic, new_whole_linker_coors, line_number)
    #     #     initial_protac_dic = protac_dic_combine(decoy_warheads_dic, warheads_atom_delete_list, new_linker_wholefile_dic, linker_atom_delete_list)
    #     #     renumber_protac_dic = protac_dic_renmuber(initial_protac_dic)
    #     #     ternary_model_dic = renumber_protac_dic.copy()
    #     #     ternary_model_dic.update(decoy_pdb_dic)

    #     #     dic_to_pdb(ternary_model_dic, ternary_model)
    #     #     ternary_model_number += 1




