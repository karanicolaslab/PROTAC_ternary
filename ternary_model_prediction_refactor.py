from itertools import product
from PROTAC_ternary import make_ternary_complex
from PROTAC_ternary import args, read_atom_file

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
    if delete_decoy_atoms != "":
        delete_decoy_atoms = [a[0] for a in read_atom_file(delete_decoy_atoms)]
    else:
        delete_decoy_atoms = []
    # label of linker overlapped atoms
    if delete_linker_atoms != "":
        delete_linker_atoms = [a[0] for a in read_atom_file(delete_linker_atoms)]
    else:
        delete_linker_atoms = []

    with open(arguments.rmsd, "w") as fwr:
        i = 0
        for d, l in product(decoys, linkers):
            rmsd, ternary_complex = make_ternary_complex(d,
                                                         l,
                                                         decoy_atoms,
                                                         linker_atoms,
                                                         delete_decoy_atoms,
                                                         delete_linker_atoms)

            if rmsd < arguments.cutoff:
                if arguments.ternary == "default":
                    out_file = "ternary{0}.pdb".format(i)
                    i += 1
                else:
                    dname = d.split("/")[-1].split(".")[0]
                    lname = l.split("/")[-1].split(".")[0]
                    out_file = dname + "_" + lname + ".pdb"
                ternary_complex.save(out_file)

            fwr.write(" ".join([d, l, str(rmsd)]) + "\n")

    fwr.close()
