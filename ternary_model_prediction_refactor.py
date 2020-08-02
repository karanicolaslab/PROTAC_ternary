from itertools import product
from PROTAC_ternary import args, read_atom_file
from PROTAC_ternary import PDBContainer


if __name__ == "__main__":
    # input arguments
    arguments = args()

    # set variables
    decoys = arguments.decoy
    linkers = arguments.linker
    decoy_atoms = arguments.decoy_alignment
    linker_atoms = arguments.linker_alignment
    delete_decoy_atoms = arguments.warheads_delete
    delete_linker_atoms = arguments.linker_delete

    # reference atoms for alignment
    decoy_atoms = read_atom_file(decoy_atoms)
    # target atoms for alignment
    linker_atoms = read_atom_file(linker_atoms)
    # label of decoys overlapped atoms
    delete_decoy_atoms = read_atom_file(delete_decoy_atoms)
    # label of linker overlapped atoms
    delete_linker_atoms = read_atom_file(delete_linker_atoms)

    fwr = open(arguments.rmsd, "w")

    i = 0

    for d, l in product(decoys, linkers):
        decoy = PDBContainer(d)
        linker = PDBContainer(l)

        decoy_atoms_ids = [decoy.atom_ids[idx] for idx in decoy_atoms]
        linker_atoms_ids = [linker.atom_ids[idx] for idx in linker_atoms]
        atoms_ids = list(zip(linker_atoms_ids, decoy_atoms_ids))

        rmsd = linker.align(decoy, atoms_ids)

        decoy.delete(delete_decoy_atoms)
        linker.delete(delete_linker_atoms)

        decoy.merge(linker)
        decoy.delete_hetHs()

        decoy_resn = [s.split(" ")[0] for s in decoy_atoms]
        linker_resn = [s.split(" ")[0] for s in linker_atoms]

        for resn in set(decoy_resn+linker_resn):
            decoy.rename(resn, "LG1", "X", 1)

        if rmsd < arguments.cutoff:
            if arguments.ternary == "default":
                out_file = "ternary{0}.pdb".format(i)
                i += 1
            else:
                dname = d.split("/")[-1].split(".")[0]
                lname = l.split("/")[-1].split(".")[0]
                out_file = dname + "_" + lname + ".pdb"

            decoy.save(out_file)

        fwr.write(" ".join([d, l, str(rmsd)]) + "\n")

    fwr.close()
