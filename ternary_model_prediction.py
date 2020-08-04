"""Wrapper code for make_ternary_complex function
"""
from itertools import product
from PROTAC_ternary import make_ternary_complex
from PROTAC_ternary import args, read_atom_file

if __name__ == "__main__":
    # get input arguments and assign to variables
    arguments = args()
    decoys = arguments.decoy
    linkers = arguments.linker
    d_alignment = arguments.decoy_alignment
    l_alignment = arguments.linker_alignment
    d_remove = arguments.warheads_delete
    l_remove = arguments.linker_delete
    rmsd_filename = arguments.rmsd
    write_complex_mode = arguments.ternary
    cutoff = arguments.cutoff

    ref_algn_atoms = read_atom_file(d_alignment)  # reference atoms of decoys
    mob_algn_atoms = read_atom_file(l_alignment)  # mobile atoms of linkers

    if d_remove != "":  # decoy atoms which will be removed after merging
        delete_decoy_atoms = [a[0] for a in read_atom_file(d_remove)]
    else:
        delete_decoy_atoms = []

    if l_remove != "":  # linker atoms which will be removed after merging
        delete_linker_atoms = [a[0] for a in read_atom_file(l_remove)]
    else:
        delete_linker_atoms = []

    with open(rmsd_filename, "w") as fwr:
        ternary_idx = 0

        # make pairs of decoy and liker files
        for d, l in product(decoys, linkers):
            # print(f"{d} and {l} have been passed to the system")

            rmsd, struct = make_ternary_complex(d, l,
                                                ref_algn_atoms,
                                                mob_algn_atoms,
                                                delete_decoy_atoms,
                                                delete_linker_atoms)

            fwr.write(" ".join([d, l, str(rmsd)]) + "\n")

            # write structure if rmsd between aligned decoy and linker
            # pass the cut off threshold
            if rmsd < cutoff:
                if write_complex_mode == "default":
                    out_file = "ternary{0}.pdb".format(ternary_idx)
                    ternary_idx += 1
                    struct.save(out_file)
                elif write_complex_mode == "specify":
                    dname = d.split("/")[-1].split(".")[0]
                    lname = l.split("/")[-1].split(".")[0]
                    out_file = dname + "_" + lname + ".pdb"
                    ternary_idx += 1
                    struct.save(out_file)

        fwr.close()

    print(str(ternary_idx) + " ternary model(s) has/have been generated.")
