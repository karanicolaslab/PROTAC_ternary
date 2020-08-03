import argparse

def args():
    parser = argparse.ArgumentParser(
                        description="""This script takes protien-protein 
                                       docking decoy(s), align with linker 
                                       conformer(s) and ouputs ternary model(s) 
                                       when alignment RMSD is less than the 
                                       cutoff""")

    parser.add_argument("-da", 
                        "--decoy_alignment",
                        required=True,
                        help="""Text file of decoy atoms for the alignment""")

    parser.add_argument("-la",
                        "--linker_alignment",
                        required=True,
                        help="""Text file of linker atoms for the alignment""")
    parser.add_argument("-d",
                        "--decoy",
                        nargs="+",
                        required=True,
                        help="""Decoys  *.pdb files""")

    parser.add_argument("-l",
                        "--linker",
                        nargs="+",
                        required=True,
                        help="""Linker conformers *.pdb files""")

    parser.add_argument("-c",
                        "--cutoff",
                        type=float,
                        default=0.4,
                        help="""Cut off value of rmsd: default is 0.4""")

    parser.add_argument("-r",
                        "--rmsd",
                        default="rmsd.txt",
                        help="""Output rmsd: rmsd values of each decoy and 
                        linker conformer with rmsd.txt as default name""")

    parser.add_argument("-t",
                        "--ternary",
                        default="default",
                        choices=["default", "specify"],
                        required=True,
                        help="""Output ternary: if not applied, no output 
                                ternary; if applied, choose either 'default' 
                                or 'specify'""")

    parser.add_argument("-wd",
                        "--warheads_delete",
                        default="",
                        help="""Warhead atoms to be deleted""")

    parser.add_argument("-ld",
                        "--linker_delete",
                        default="",
                        help="""Linker atoms to be deleted""")

    parser.add_argument("-v",
                        "--version",
                        action="version",
                        version="Version 1.0")

    args = parser.parse_args()

    if args.cutoff < 0.0:
        print("Error: cut_off must be positive.")

    return args


def read_atom_file(file):
    atoms = [line.rstrip('\n').split() for line in open(file,"r")]
    atoms = [[a[0] + " " + aa for aa in a[1::]] for a in atoms ]

    return atoms