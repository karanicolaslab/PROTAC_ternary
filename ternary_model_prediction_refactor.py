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
    delete_decoy_atoms = [a[0] for a in read_atom_file(delete_decoy_atoms)]
    # label of linker overlapped atoms
    delete_linker_atoms = [a[0] for a in read_atom_file(delete_linker_atoms)]

    fwr = open(arguments.rmsd, "w")

    i = 0

    for d, l in product(decoys, linkers):
        decoy = PDBContainer(d)
        linker = PDBContainer(l)

        rmsd = linker.align(decoy, targ_atoms=linker_atoms, ref_atoms=decoy_atoms)

        decoy.delete_atoms(delete_decoy_atoms)
        linker.delete_atoms(delete_linker_atoms)

        decoy.merge(linker)
        decoy.delete_hetHs()

        decoy_resn = [ss.split(" ")[0] for s in decoy_atoms for ss in s]
        linker_resn = [ss.split(" ")[0] for s in linker_atoms for ss in s]

        for resn in set(decoy_resn + linker_resn):
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
