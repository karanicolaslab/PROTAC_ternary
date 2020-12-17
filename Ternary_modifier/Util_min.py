import argparse


def args():

    parser = argparse.ArgumentParser()
    parser.add_argument("-s",
                        "--input_pdb_list",
                        required=True)

    args = parser.parse_args()

    return args


                        
