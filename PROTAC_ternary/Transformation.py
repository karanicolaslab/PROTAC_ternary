import numpy as np
import math

def apply_rotation(rot1, rot2, rot3, orig_coors):    
    new_coors = np.copy(orig_coors)
    # consolidate the three rotations into a single transformation matrix
    rot_mat1 = np.array([[1, 0, 0],
                         [0, math.cos(rot1), -math.sin(rot1)],
                         [0, math.sin(rot1), math.cos(rot1)]])

    rot_mat2 = np.array([[math.cos(rot2), 0, math.sin(rot2)],
                         [0, 1, 0],
                         [-math.sin(rot2), 0, math.cos(rot2)]])

    rot_mat3 = np.array([[math.cos(rot3), -math.sin(rot3), 0],
                         [math.sin(rot3),  math.cos(rot3), 0],
                         [0, 0, 1]])

    rot_mat_all = np.dot(np.dot(rot_mat1, rot_mat2), rot_mat3)

    # apply the rotation matrix
    for index in range(len(new_coors[:, 0])):
        new_coors[index, :] = rot_mat_all.dot(new_coors[index, :])

    return new_coors


def apply_tranlation(xtrans, ytrans, ztrans, orig_coors):
    """
       Take in original coordinates and transformation vetors, and output new 
       coordinates after transformation. Transformation will consist of first 
       rotating, then translating
    """

    new_coors = np.copy(orig_coors)
    new_coors[:, 0] += xtrans
    new_coors[:, 1] += ytrans
    new_coors[:, 2] += ztrans

    return new_coors

# def apply_transformation(rot1, rot2, rot3, xtrans, ytrans, ztrans, orig_coords):
#     rot_orig_coords = apply_rotation(rot1, rot2, rot3, orig_coords)
#     trans_orig_coords = apply_tranlation(xtrans, ytrans, ztrans, rot_orig_coords)

#     return trans_orig_coords


def apply_transformation(rot1, rot2, rot3, xtrans, ytrans, ztrans, orig_coords):
    """Take in original coordinates and transformation vetors, and output new coordinates after transformation."""
    # transformation will consist of first rotating, then translating
    new_coors = np.copy(orig_coords)
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