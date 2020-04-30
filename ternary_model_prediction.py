import math
import numpy as np
import collections as col
import scipy.optimize
import argparse


# note: open list file to get decoy/conformer list
def read_file(name_file):
	name_list = []
	open_file = open(name_file, "r")
	for line in open_file:
		name = line.strip()
		name_list.append(name)
	open_file.close()
	return name_list

# note: atom files to list and the first element is the residue ID, for both alignment atoms and delete atoms.
def atom_list_generation(atom_file):
	atom_list = [line.rstrip('\n').split() for line in open(atom_file)]
	return atom_list

# note: transfer linker/decoy pdb to a dictionary with (ligand name + atom name) as key and (x, y, z) as val 
def pdb_to_dictionary(pdb_file):
	pdb_dic = col.OrderedDict()
	f_open = open(pdb_file, "r")
	for line in f_open:
		if line.startswith("HETATM"):
			residue_name = line[17:20].strip()
			atom_name = line[12:16].strip()
			key = residue_name + " " + atom_name
			val = line
			pdb_dic[key] = val
			#i += 1
		if line.startswith("ATOM"):
			key = line[0:26]
			val = line
			pdb_dic[key] = val
	return pdb_dic

# note: creat coordinates array (6X3)
def pdb_coors_assignment(pdb_list, pdb_dic):
	pdb_coors = np.empty([6,3])
	for i in range(len(pdb_list)):
		atom = []
		for j in range(len(pdb_list[i])):
			atom.append(pdb_list[i][j])
		x_coor = []
		y_coor = []
		z_coor = []
		for n in range(len(atom)-1):
			m = n + 1
			ligand_name = atom[0]
			atom_name_list = atom[m]
			key = ligand_name + " " + atom_name_list
			coors_line = pdb_dic[key]
			x_coor.append(float(coors_line[30:38].strip()))
			y_coor.append(float(coors_line[38:46].strip()))
			z_coor.append(float(coors_line[46:54].strip()))

		pdb_coors[i,0] = np.average(x_coor)
		pdb_coors[i,1] = np.average(y_coor)
		pdb_coors[i,2] = np.average(z_coor)
	return pdb_coors

def count_line (linker_dic):
	line_number = 0
	for key in linker_dic:
		line_number += 1
	return line_number

# note: assigne whole linker coors
def whole_linker_coors_assignment(linker_dic, atom_number):
	whole_linker_conformer_coors = np.empty([atom_number,3])
	n = 0
	for key in linker_dic:
		value = linker_dic[key]
		if n < atom_number:
			whole_linker_conformer_coors[n,0] = float(value[30:38])
			whole_linker_conformer_coors[n,1] = float(value[38:46])
			whole_linker_conformer_coors[n,2] = float(value[46:54])
			n += 1
		else:
			break
	return whole_linker_conformer_coors

def apply_transformation (xtrans_ytrans_ztrans_rot1_rot2_rot3, orig_coors):
    # transformation will consist of first rotating, then translating
	new_coors = np.copy(orig_coors)
	xtrans, ytrans, ztrans, rot1, rot2, rot3 = xtrans_ytrans_ztrans_rot1_rot2_rot3
    # consolidate the three rotations into a single transformation matrix
	rot_mat1 = np.array([[ 1, 0, 0], [ 0, math.cos(rot1), -math.sin(rot1)], [ 0, math.sin(rot1), math.cos(rot1)]])
	rot_mat2 = np.array([[ math.cos(rot2), 0, math.sin(rot2)], [ 0, 1, 0], [ -math.sin(rot2), 0, math.cos(rot2)]])
	rot_mat3 = np.array([[ math.cos(rot3), -math.sin(rot3), 0], [ math.sin(rot3), math.cos(rot3), 0], [ 0, 0, 1]])
	rot_mat_all = np.dot(np.dot(rot_mat1, rot_mat2), rot_mat3)

    # apply the rotation matrix
	for index in range(len(new_coors[:,0])):
		new_coors[index,:] = rot_mat_all.dot(new_coors[index,:])

    # apply the translation
	new_coors[:,0] += xtrans
	new_coors[:,1] += ytrans
	new_coors[:,2] += ztrans
	return new_coors

    
# note: this is not actually the RMSD, it's the sum of the squares of the error (for speed)
# note: to get rmsd, divide by number of atoms and then take the square root
def eval_rmsd_after_transformation(xtrans_ytrans_ztrans_rot1_rot2_rot3, moving_coors, ref_coors):
	xtrans, ytrans, ztrans, rot1, rot2, rot3 = xtrans_ytrans_ztrans_rot1_rot2_rot3
	new_coors = apply_transformation(xtrans_ytrans_ztrans_rot1_rot2_rot3, moving_coors)
	dist = new_coors - ref_coors
	dist = dist * dist
	sq_err = np.sum(dist)
	return sq_err

# note: store the lines of the new conformer after linker_rmsd optimization
def linker_conformer_after_rmsd_optimizaton(initial_linker_dic, new_wholelinker_coors_array, linker_atom_number):
	new_linker_dic = col.OrderedDict()
	m = 0
	for key in initial_linker_dic:
		value = initial_linker_dic[key]
		newkey = key
		if m < linker_atom_number:
			x = str('%8.3f' % (float(new_wholelinker_coors_array[m,0]))).rjust(8)
			y = str('%8.3f' % (float(new_wholelinker_coors_array[m,1]))).rjust(8)
			z = str('%8.3f' % (float(new_wholelinker_coors_array[m,2]))).rjust(8)
			newval_x = value.replace(value[30:38], x)
			newval_y = newval_x.replace(newval_x[38:46], y)
			newval_z = newval_y.replace(newval_y[46:54], z)
			new_linker_dic[newkey] = newval_z
			m += 1
		else:
			new_linker_dic[newkey] = value
	return new_linker_dic

def dic_split(input_dic, split_keyword):
	new_dic = col.OrderedDict()
	for key in input_dic:
		if input_dic[key].startswith(split_keyword):
			val = input_dic[key]
			new_dic[key] = val
	return new_dic

def protac_dic_combine (warheads_dic, warheads_atom_delete_list, linker_dic, linker_atom_delete_list):
	protac_dic = col.OrderedDict()
	warheads_dic_inter = warheads_dic.copy()
	linker_dic_inter = linker_dic.copy()
	for item in warheads_atom_delete_list:
		key = item[0] + " " + item[1]
		del warheads_dic_inter[key]
	for item in linker_atom_delete_list:
		key = item[0] + " " + item[1]
		del linker_dic_inter[key]
	warheads_dic_inter.update(linker_dic_inter)
	protac_dic = warheads_dic_inter.copy()
	for key in warheads_dic_inter:
		check = (str(key).split()[1]).strip()
		if check.startswith("H"):
			del protac_dic[key]
	#protac_dic = warheads_dic_inter.copy()
	return protac_dic

def protac_dic_renmuber (old_protac_dic):
	new_protac_dic = col.OrderedDict()
	i = 1
	for key in old_protac_dic:
		atom_number = str(i).rjust(5)
		line = old_protac_dic[key]
		newline = line[0:6] + atom_number + line[11:14] + "   LGX X" + line[22:]
		new_protac_dic[key] = newline
		i += 1
	return new_protac_dic

def dic_to_pdb (dic, pdb_name):
	f_open = open(pdb_name, "w")
	for key in dic:
		line = str(dic[key])
		#print line
		f_open.write(line)
	f_open.close()
	return 0


def main():

	# note: decoy & linker would be store in the list file
	decoy_list =[]
	linker_list = []

	# reading flags and input files
	parser = argparse.ArgumentParser(
		description='This script takes protien-protein docking decoy(s), align with linker conformer(s) and ouputs ternary model(s) when alignment RMSD is less than the cutoff.')
	parser.add_argument("-da", "--decoy_alignment", 
		help = "Text file of decoy atoms for the alignment")
	parser.add_argument("-la", "--linker_alignment", 
		help = "Text file of linker atoms for the alignment")
	parser.add_argument("-d", "--decoy", 
		help = "Decoys, either one *.pdb or a file with multiple lines: /path/to/decoy.pdb")
	parser.add_argument("-l", "--linker", 
		help = "Linker conformers, either one *.pdb or a file with multiple lines: /path/to/linker.pdb")
	parser.add_argument("-c", "--cutoff", type = float, default = 0.4, 
		help = "Cut off value of rmsd: default is 0.4")
	parser.add_argument("-r", "--rmsd", default = "rmsd.txt", 
		help = "Output rmsd: rmsd values of each decoy and linker conformer with rmsd.txt as default name")
	parser.add_argument("-t", "--ternary", default = "",  
		help = "Output ternary: if not applied, no output ternary; if applied, choose either 'default' or 'specify'")
	parser.add_argument("-wd", "--warheads_delete", default = "",
		help = "Warhead atoms to be deleted")
	parser.add_argument("-ld", "--linker_delete", default = "",
		help = "Linker atoms to be deleted")
	parser.add_argument("-v", "--version", 
		help = "Version 1.0", action="store_true")
	
	args = parser.parse_args()

	# check if input files are correct
	if args.version:
		print("\nThis is version 1.0 of the ternary_model_prediction code.\n")
		quit()
	if not args.decoy_alignment:
		print ("Error: must specify a file containing the list of decoy atoms which would be used for alignment.")
		quit()
	if not args.linker_alignment:
		print ("Error: must specify a file containing the list of linker atoms which would be used for alignment.")
		quit()
	if not args.decoy:
		print("Error: must specify an input decoy pdb / decoy list file.")
		quiy()
	if not args.linker:
		print("Error: must specify an input linker pdb / linker list file.")
		quit()
	if not args.ternary or args.ternary == "":
		print("No ternary structure would be genrated.")
	if not args.cutoff or args.cutoff == 0.4:
		print("No customized cut off value, default value(0.4) will be applied.")
	if not args.rmsd or args.rmsd == "rmsd.txt":
		print("No customized rmsd file name, default name(rmsd.txt) will be applied.")
	if not args.warheads_delete or args.warheads_delete == "":
		print("Warning: all atoms from warheads will be included in the protac.")
	if not args.linker_delete:
		print("Warning: all atoms from linker conformer will be included in the protac.")


	linker_atoms_alignment = args.linker_alignment
	decoy_atoms_alignment = args.decoy_alignment

	# input decoy: pdb or list file
	decoy_file = args.decoy
	if (decoy_file[-4:] == ".pdb"):
		decoy_list.append(decoy_file)
	else:
		decoy_list = read_file(decoy_file)

	# input linker: pdb or list file
	linker_file = args.linker
	if (linker_file[-4:] == ".pdb"):
		linker_list.append(linker_file)
	else:
		linker_list = read_file(linker_file)

	# ternary structure: generate or not; if yes, what kind of name
	ternary_model_name = args.ternary
	if (ternary_model_name == "default"):
		print ("Default name will be applied to ternary model.")
	elif (ternary_model_name == "specify"):
		print ("Specified name (decoy + linker) will be applied to ternary model.")
	elif ternary_model_name == "":
		pass
	else:	
		print ("Error: unrecognized setting of '-t/--ternary':", ternary_model_name, "\n-t/--ternary can be either 'default' or 'specify'.")
		quit()


	cut_off = args.cutoff
	if (cut_off < 0):
		print ("Error: cut_off must be positive.")
		quit()

	rmsd_file = args.rmsd

	warheads_atom_delete = args.warheads_delete
	if (warheads_atom_delete == ""):
		warheads_atom_delete_list = []
	else:
		warheads_atom_delete_list = atom_list_generation(warheads_atom_delete)

	linker_atom_delete = args.linker_delete
	if (linker_atom_delete == ""):
		linker_atom_delete_list = []
	else:
		linker_atom_delete_list = atom_list_generation(linker_atom_delete)

	

	decoy_atom_alignment_list = atom_list_generation(decoy_atoms_alignment)
	linker_atom_alignment_list = atom_list_generation(linker_atoms_alignment)

	decoy_array = []
	decoy_file_array = []
	for decoy_name in decoy_list:
		decoy_dic = pdb_to_dictionary(decoy_name)
		#print len(decoy_dic)
		file_array_item = [decoy_name, decoy_dic]
		decoy_file_array.append(file_array_item)
		docking_decoy_coors = pdb_coors_assignment(decoy_atom_alignment_list, decoy_dic)
		array_item = [decoy_name, docking_decoy_coors]
		decoy_array.append(array_item)
	#print decoy_array[0][1][0]
	print (str(len(decoy_list)) + " decoy(s) has/have been passed to the system.")

	linker_array = []
	linker_file_array = []
	for linker_name in linker_list:
		linker_dic = pdb_to_dictionary(linker_name)
		file_array_item = [linker_name, linker_dic]
		linker_file_array.append(file_array_item)
		orig_linker_coors = pdb_coors_assignment(linker_atom_alignment_list, linker_dic)
		array_item = [linker_name, orig_linker_coors]
		linker_array.append(array_item)
	print (str(len(linker_list)) + " linker conformer(s) has/have been passed to the system.")


	list_to_write = []
	ternary_model_number = 0
	for i in range (len(decoy_array)):
		decoy_pdb = decoy_array[i][0]
		decoy_wholefile_dic = decoy_file_array[i][1]
		decoy_pdb_dic = dic_split(decoy_wholefile_dic, "ATOM")
		decoy_warheads_dic = dic_split(decoy_wholefile_dic, "HETATM")
		docking_decoy_coors = decoy_array[i][1]
		for j in range (len(linker_array)):
			linker_conformer_pdb = linker_array[j][0]
			linker_wholefile_dic = linker_file_array[j][1]
			#print linker_wholefile_dic
			orig_linker_coors = linker_array[j][1]
			natoms = len(orig_linker_coors[:,0])
			res = scipy.optimize.minimize(eval_rmsd_after_transformation, (0., 0., 0., 0., 0., 0.), args=(orig_linker_coors, docking_decoy_coors), method='Powell')
			new_linker_coors = apply_transformation((res.x[0], res.x[1], res.x[2], res.x[3], res.x[4], res.x[5]), orig_linker_coors)
			rmsd_value = math.sqrt(eval_rmsd_after_transformation((0., 0., 0., 0., 0., 0.), new_linker_coors, docking_decoy_coors)/natoms)
			rmsd_line = "Final RMSD between " + decoy_pdb + " and " + linker_conformer_pdb + " is " + str(rmsd_value) + "\n"
			#print rmsd_line
			write_line = decoy_pdb + " " + linker_conformer_pdb + " " + str(rmsd_value) + "\n"
			list_to_write.append(write_line)


			if (rmsd_value > cut_off or ternary_model_name == ""):
				pass
			else:
				if (ternary_model_name == "default"):
					ternary_model = "ternary" + str(ternary_model_number) + ".pdb"
				else:
					ternary_model = (decoy_pdb.split("/")[-1]).split(".")[0] + "_" + linker_conformer_pdb.split("/")[-1]

				line_number = count_line(linker_wholefile_dic)
				#print line_number 
				orig_whole_linker_coors = whole_linker_coors_assignment(linker_wholefile_dic, line_number)
				new_whole_linker_coors = apply_transformation((res.x[0], res.x[1], res.x[2], res.x[3], res.x[4], res.x[5]), orig_whole_linker_coors)
				new_linker_wholefile_dic = linker_conformer_after_rmsd_optimizaton(linker_wholefile_dic, new_whole_linker_coors, line_number)
				#print new_linker_wholefile_dic
				#print warheads_atom_delete_list, linker_atom_delete_list
				initial_protac_dic = protac_dic_combine(decoy_warheads_dic, warheads_atom_delete_list, new_linker_wholefile_dic, linker_atom_delete_list)
				#print initial_protac_dic
				renumber_protac_dic = protac_dic_renmuber(initial_protac_dic)
				ternary_model_dic = renumber_protac_dic.copy()
				#print len(ternary_model_dic)
				#print len(decoy_wholefile_dic)
				#print len(decoy_pdb_dic)
				ternary_model_dic.update(decoy_pdb_dic)
				#print len(ternary_model_dic)
				
				dic_to_pdb(ternary_model_dic, ternary_model)
				ternary_model_number += 1

				#print renumber_protac_dic
				#dic_to_pdb(initial_protac_dic, "initial_protac_test_test.pdb")
				#dic_to_pdb(renumber_protac_dic, "renumber_protac_test.pdb")

		#break
	if (rmsd_file == "none"):
		list_to_write = []
	else:
		f_rmsd = open(rmsd_file, "w")
		for item in list_to_write:
			f_rmsd.write(item)
		f_rmsd.close()

	print (str(ternary_model_number) + " ternary model(s) has/have been generated.")

main()
