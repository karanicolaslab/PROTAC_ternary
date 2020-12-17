#script to minimize full ternary complex

import os, sys
import collections
import time
import subprocess, os
from Util_min import args

BABEL = os.getenv('BABEL')
MOL2PARAMS = os.getenv('MOL2PARAMS')

if BABEL is None:
   print("Need babel path")
   sys.exit()

if MOL2PARAMS is None:
   print("Need Rosetta mol2params path")
   sys.exit()

arguments = args()


#first, read in ternary.pdb list 
pdb_lst = []  
pdb_file = open(arguments.input_pdb_list, 'r')
for line in pdb_file:
    name = line.strip()
    pdb_lst.append(name)
pdb_file.close


#next, for each terary.pdb, copy the ligand coordinates (dictionary?) into a new .pdb file
def lig_to_dictionary(ter_file):
    lig_dic = collections.OrderedDict()
    f_file = open(ter_file , "r")
    counter = 0
    for line in f_file:
        if line.startswith("HETATM"):
#            residue_name = line[17:20].strip()
#            atom_name = line[10:15].strip()
#            key = residue_name + " " + atom_name #couldn't I just do line number here instead? 
            key = counter
            val = line.rstrip('\n')
            lig_dic[key] = val
            #i+=1
            counter += 1 
        if line.startswith("ATOM"):
            pass 
    return lig_dic 

def lig_pdb_maker(ter_file,lig_dic):
    name = ter_file.strip(".pdb") + "_lig.pdb"
    f = open(name, "w")
    for value in lig_dic.values():
        f.write(value)
        f.write(os.linesep)
    f.close()
    print(name + " was succesfully created!")
    return f 

def babel_runner(ter_file):
    name = ter_file.strip(".pdb") + "_lig.pdb"
    mol2_name = name.strip(".pdb") + ".mol2"
    command = BABEL + " -h -ipdb " + name + " -omol2 -O " + mol2_name 
    print("About to run: " + command)
    os.system(command)
    return 0 

def mol2_to_params(ter_file):
    name = ter_file.strip(".pdb") + "_lig.mol2"
    command = "python2.6 " + MOL2PARAMS + " " + name + " -n TRN --clobber --no-param -p "+ ter_file.split(".")[0]
    print("About to run: " + command)
    os.system(command)
    return 0 


def insert_lig(pdb2_lig, ter_file):
    counter1 = 1
    counter2 = 1
    dic_file = pdb2_lig 
    lig2_dic = collections.OrderedDict()
    for line in dic_file:
        key = counter1
        val = line 
        lig2_dic[key] = val
        counter1 += 1
    counter1 -= 1
    lig2_dic.pop(int(counter1))
    f = open(ter_file, "r")
    mod_pdb = ter_file.strip(".pdb") + "_mod.pdb"
    ff = open(mod_pdb, 'w')
    judgement = False 
    for line in f:
        if line.startswith("ATOM"):
            judgement = True
        if line.startswith("HETATM"):
            judgement =  False
        if line.startswith("TER"):
            judgement = False 
        elif judgement:
            ff.write(line)
    f.close()
    ff.close()
    with open(mod_pdb, 'r+') as fff:
        for key in lig2_dic:
            new_line = str(lig2_dic[key])
            fff.write(new_line)
    fff.close()
    #Add in a line here where it goes back and deletes all the _mod files so it is a bit cleaner, oops.
    print(ter_file + " was modified?")
    return fff

    
def main():
    for pdb in pdb_lst:
        orig_lig  = lig_to_dictionary(pdb)
        new_lig = lig_pdb_maker(pdb, orig_lig)
    print("Moving on to generating params")
    for pdb in pdb_lst:
        mol2_lig = babel_runner(pdb)
        mol2_to_params(pdb)
        lig_file = open(pdb.split(".")[0]+"_0001.pdb", 'r')
        fin_pdb = insert_lig(lig_file, pdb)


main()
    

#Just extra commands that I didn't need but could come in handy later in time
#    p = subprocess.Popen(['/bin/bash', '-c', command])                                                                                                              
# p.wait()
#            cmd = "sed -e " + str(counter2) + "d " + ter_file +        " >> " + ter_file.strip(".pdb") + "_mod" + str(counter2) + ".pdb"
#            os.system(cmd)
#            pass
#        if line.startswith("HETATM"):
#            cmd2 = "sed -e " + str(counter2) + "d " + ter_file.strip(".pdb") + "_mod" + str(counter2) + ".pdb" + " >> " + ter_file.strip(".pdb") + "_mod" + str(counter2) + ".pdb"
#            os.system(cmd2)
#            counter2 += 1
#            continue
#        else:
#            pass~
