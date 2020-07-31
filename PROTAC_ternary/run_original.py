cd test/
python ../ternary_model_prediction.py -d CDK6_CRBN_docked_pdbs.txt -l YKL_06_102_conformers.txt -da decoy_atom_list.txt -la linker_atom_list.txt -t default -wd decoy_atom_delete.txt -ld linker_atom_delete.txt
cd ../