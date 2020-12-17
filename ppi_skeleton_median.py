### Isolates the skeleton I_sc score and calculates the median value from minimization output files

import os
import pandas as pd

def ppi_scores_file(out_file):
    with open(out_file, 'r') as output_file:
        with open("ppi_scores.txt", 'w') as scores:
            scores.write("Interface_Scores:TAG  input_pdb_name  bound_energy Interface_Energy Total_BSA Interface_HB  Total_packstats Interface_unsat \n")
            for line in output_file:
                if line.startswith("Interface_Scores:") and "TAG" not in line:
                    scores.write(line)

def median_calculator():
    sc = pd.read_csv("ppi_scores.txt", delim_whitespace=True )
    ff = open("skeleton_median.txt", 'w')
    print(sc.median(axis = 0))
    output = str(sc.median(axis = 0))
    ff.write(output)
    ff.close()
                
def main():
    out_file = input(str("What is the name of the output file? "))
    ppi_scores_file(out_file)
    median_calculator()
            





main()
