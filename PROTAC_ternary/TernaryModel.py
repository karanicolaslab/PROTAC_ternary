from PROTAC_ternary import PDBContainer


def make_ternary_complex(decoy_file,
                         linker_file,
                         decoy_atoms,
                         linker_atoms,
                         aln_iters,
                         delete_decoy_atoms=[],
                         delete_linker_atoms=[]):
    """The main function for building ternary complex

       At the first step the function will align linker position relative
       warhead positions and then merge both structure to one. If
       delete_decoy_atoms and delete_linker_atoms are presented, the function
       will remove these atoms as well as all HET hydrogens atoms.

    Args:
        decoy_file (str): File name of warhead PDB
        linker_file (str): File name of linker PDB
        decoy_atoms (list): List of warheads atoms that will be used for
                            alignment
        linker_atoms (list): List of linker atoms that will be used for
                             alignment
        delete_decoy_atoms (list, optional): list of warheads atoms that will
                                             be removed after merging
        delete_linker_atoms (list, optional): list of linkers atoms that will
                                              be removed after merging

    Returns:
        (float, PDBContainer): aligned and merged protein/
                               /warheads/linker structure
    """
    decoy = PDBContainer(decoy_file)
    linker = PDBContainer(linker_file)

    rmsd = linker.align(decoy, 
                        targ_atoms=linker_atoms,
                        ref_atoms=decoy_atoms,
                        aln_iters=aln_iters)

    if len(delete_decoy_atoms) != 0:
        decoy.delete_atoms(delete_decoy_atoms)

    if len(delete_linker_atoms) != 0:
        linker.delete_atoms(delete_linker_atoms)

    decoy.merge(linker)
#    decoy.delete_hetHs()

    decoy_resn = [ss.split(" ")[0] for s in decoy_atoms for ss in s]
    linker_resn = [ss.split(" ")[0] for s in linker_atoms for ss in s]

    for resn in set(decoy_resn + linker_resn):
        decoy.rename(resn, "LG1", "X", 1)

    return rmsd, decoy
