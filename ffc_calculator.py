#Program to calculate the FFC value based on user input
#NOTE: the ternary models must all have their minimization scores put into a single score.sc file for this to work

import os
import pandas as pd


def filtered_complexes(sc_file):#add file as argument
    skel_med = float(input("Enter the median interaction energy for the docked cohort without the PROTAC linker: "))
    sc = pd.read_csv(sc_file, delim_whitespace=True)
    df = sc[sc.Interface_Energy < skel_med]
    fc = len(df.index)
    print(df)
    return df[['input_pdb_name', 'Interface_Energy']]    
 
def ffc_eqn(dm, lcm, fc):
    ffc = (fc*1000)/(dm*lcm)
    return ffc
    
def main():
    out_file = input(str("What is the name of the ternary csv output file?"))
    dm = int(input("Docked models created: "))
    lcm = int(input("Linker conformers created: "))
    fc = filtered_complexes(out_file)
    ffc = ffc_eqn(dm, lcm, fc.shape[0])
    print(fc.shape[0],"passed the energy filter.")
    fc.to_csv("enery_passed.txt",index=False)
    print("FFC value is ", ffc)

main()
