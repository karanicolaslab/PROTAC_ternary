"""This script takes in the ternary models created from the ternary_model_prediction.py and the .params file created from the PROTAC of the model with the highest RMSD value. It then uses the Rosetta ppi_minimize function to minimize all the models with this .params file to output interaction scores via the output file. These scores will be used for the energy filter and then to subsequently calculate the FFC value for this given ternary complex.
"""
import os

def pdb_names(file):
    """Reads in the list of ternary model pdb files

    Args:
        file (list) : list of pdbs

    Returns:
        (list) : list of pdbs stripped
    """
    with open(file, 'r') as pdb_file:
        pdb_lst = [line.strip() for line in pdb_file]

    return pdb_lst


def get_hetatoms(ter_file, onlyHet=True):
    """Extracts the PROTAC from the ternary model

    Args:
        ter_file (file): TernaryModel PDB file
        onlyHet (bool, optional): Atoms from the PROTAC ligand

    Returns:
        (list) : A list of all the PROTAC's atoms
    """
    atoms = []
    with open(ter_file, "r") as f_file:
        for line in f_file:
            atoms.append(line.rstrip('\n'))

    if onlyHet:
        return []

    if onlyHet:
        atoms = filter(lambda x: x.startswith("HETATM"), atoms)

#    if onlyAtm:
#        atoms = filter(lambda x: x.startswith("ATOM"), atoms)

    return list(atoms)


def pdb_maker(name, atoms):
    """Writes the pdb file for the PROTAC 

    Args:
        name (str): The name of the ternary model ligand file (PROTAC)
        atoms (list): The list of the PROTAC's atoms
    """
    with open(name, "w") as f:
        text = "\n".join(atoms)
        f.write(text)

    f.close()


def babel_runner(pdb_file, mol_file):
    """Runs obabel to generate a mol2 file of the PROTAC 

    Args:
        pdb_file (str): name of the PROTAC-only pdb file
        mol_file (str): name of the PROTAC lig mol2 file
    """
#    BABEL = "/mnt/shared_applications/openbabel/2.4.1/bin/obabel"
    cmd = f"obabel -h -ipdb {pdb_file} -omol2 -O {mol_file}"
    os.system(cmd)


def mol2_to_params(ter_file, lig_name="TRN"):
    """Use the Rosetta script to convert the mol2 file to a pdb file compatible with the PARAMS file in minimization

    Args:
        ter_file (str): name of the PROTAC lig mol2 file 
    """
    # do we want the user to specify the three letter identifier to match the
    # inputed params value earlier? I think we should!
    global MOL2PARAMS
    cmd = f"{MOL2PARAMS} {ter_file} -n {lig_name} --clobber --no-param"
    os.system(cmd)


def insert_lig(ter_file, mod_pdb):
    """Take the pdb file from the mol2_to_params function and place it inside the original ternary model

    Args:
        ter_file (str) : name of original ternary complex pdb file
        mod_pdb (str) : name of modified ternary complex pdb file
    """
    protein_atoms = get_atoms(ter_file, onlyAtm=True)

    complex_atoms = protein_atoms + ["TRN_0001.pdb"]

    with open(mod_pdb, 'w') as f:
        f.write("\n".join(complex_atoms))
        f.close()


def rosetta_minimizer(params, ter_file):
    """Summary

    Args:
        params (file): Params file taken from the pose with the lowest rmsd 
        value ter_file (str): name of the modified pdb file with updated PROTAC
    """
    global MINIMIZER_APP
    cmd = f"{MINIMIZER_APP} -s {mod_pdb} -extra_res_fa" + f" {params} -jump_all"
    os.system(cmd)


def remove_file(file):
    """Remove the generated PROTAC ligand file from the directory

    Args:
        file (str): name of the TRN_0001.pdb that was generated from
        mol_2_params function
    """
    cmd = f"rm {file}"
    os.system(cmd)


def main():
    """Main function that inputs the list of ternary models and the best rmsd
       params file and performs minimizations.
    """
    pdb_lst = pdb_names("pdb_list.txt")
    params = "TRN.params"

    for pdb in pdb_lst:
        hetatms = get_hetatoms(pdb, onlyHet=True)  # extract hetatoms from PDB

        lig_filename = pdb.strip(".pdb") + "_lig.pdb"
        pdb_maker(lig_filename, hetatms)  # save hetatoms to other file

        mol2_lig_filename = lig_filename.strip(".pdb") + ".mol2"
        babel_runner(lig_filename, mol2_lig_filename)  # convert to mol2
        mol2_to_params(mol2_lig_filename)  # generate new pdb 

        # skip iteration if mol2params was failed
        if not os.path.exists(str(os.getcwd()) + "/TRN_0001.pdb"):
            continue

        mod_pdb = pdb.strip(".pdb") + "_mod.pdb"
        insert_lig(pdb, mod_pdb)  # insert ligand into protein

        remove_file(lig_filename)  # remove ligand pdb

        mod_pdb = ter_file.strip(".pdb") + "_mod.pdb"
        rosetta_minimizer(params, mod_pdb)  # minimize ternary complex

        remove_file("TRN_0001.pdb")  # remove ligand pdb


if __name__ == "__main__":
    
    if "ROSETTA_DIR" not in os.environ:
        raise ValueError("Please, add ROSETTA_DIR into environment")
    else:
        ROSETTA_DIR = os.environ.get('ROSETTA_DIR')
        MINIMIZER_APP = f"{ROSETTA_DIR}/main/source/bin/minimize_ppi.default.linuxgccrelease"
        MOL2PARAMS = f"{ROSETTA_DIR}/main/source/scripts/python/public/molfile_to_params_yates.py"
    
    main()
