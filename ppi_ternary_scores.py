###Isolates I_sc scores from minimiztion output files to be used as input for ffc calculation input                                                      
import os

def ppi_scores_file(out_file):
    with open(out_file, 'r') as output_file:
        with open("ppi_ternary_scores.csv", 'w') as scores:
            scores.write("Interface_Scores:TAG  input_pdb_name  bound_energy Interface_Energy Total_BSA Interface_HB  Total_packstats Interface_unsat \n")
            for line in output_file:
                if line.startswith("Interface_Scores:") and "TAG" not in line:
                    scores.write(line)


def main():
    out_file = input(str("What is the name of the output file? "))
    ppi_scores_file(out_file)


main()
