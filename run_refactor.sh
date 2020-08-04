# Test general functionality
echo "Test #1"
python ternary_model_prediction.py -d test/CDK6_CRBN_dock*pdb \
                                   -l test/YKL_06_102_conformer*pdb \
                                   -da test/decoy_atom_list.txt \
                                   -la test/linker_atom_list.txt \
                                   -c 10.0 \
                                   -r rmsd.txt \
                                   -t specify \
                                   -wd test/decoy_atom_delete.txt \
                                   -ld test/linker_atom_delete.txt
   