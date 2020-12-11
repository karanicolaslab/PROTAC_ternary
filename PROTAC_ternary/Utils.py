import argparse


def args():
    """Process input arguments

    Returns:
        argparse.Namespace: Namespace of arguments
    """
    parser = argparse.ArgumentParser(
        description="""This script takes protien-protein docking decoy(s),
                       align with linker conformer(s) and ouputs ternary
                       model(s) when alignment RMSD is less than the cutoff""")

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
                        help="""Decoys, one or several PDB files""")
    parser.add_argument("-dl",
			"--decoy_list",
			help="""Decoy PDB files listed in a single txt file in PWD""") 

    parser.add_argument("-l",
                        "--linker",
                        nargs="+",
                        help="""Linker, one or several PDB files""")

    parser.add_argument("-ll",
			"--linker_list",
			help ="""Linker PDB files listed ina signle text file in PWD""")

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
                        default="none",
                        choices=["none", "default", "specify"],
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

    parser.add_argument("-ai",
                        "--alignment_iterations",
                        type=int,
                        default=5,
                        help="""Number of iterations used for alignment""")

    parser.add_argument("-v",
                        "--version",
                        action="version",
                        version="Version 1.0")

    args = parser.parse_args()
    print (args)

    if args.cutoff < 0.0:
        print("Error: cut_off must be positive.")
    if args.linker is not None and args.linker_list is not None:
        raise Exception ("Linker and Linker List cannot both be defined")
    elif args.linker is None and args.linker_list is None:
        raise Exception ("Linker and Linker_List cannot both be empty")
    if args.decoy is not None and args.decoy_list is not None:
        raise Exception ("Decoy and Decoy List cannot both be defined")	
    elif args.decoy is None and args.decoy_list is None: 
        raise Exception ("Decoy and Decoy_List cannot both be empty")  
    if args.ternary == "default":
        print("Default name will be applied to ternary model")
    elif args.ternary == "specify":
        print("Specified name (decoy+linker) will be applied to ternary model")
    if args.decoy_list is not None:
        args.decoy=read_list_file(args.decoy_list)
    if args.linker_list is not None:
        args.linker=read_list_file(args.linker_list)
    print (args)
    return args
   

def read_list_file(list_file):
    with open (list_file,"r") as file:
        file=[i for i in file.read().split("\n") if i!= "" ]
        return file 
def read_atom_file(file):
    """Read file with atoms label

    Args:
        file (str): Name of file with atoms labels

    Returns:
        list: atoms labeled as <Residue Name> <Atom Name>
    """
    atoms = [line.rstrip('\n').split() for line in open(file, "r")]
    atoms = [[a[0] + " " + aa for aa in a[1::]] for a in atoms]

    return atoms
