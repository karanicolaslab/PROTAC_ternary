# # Test general functionality
# echo "Test #1"
# python ternary_model_prediction_refactor.py -d test/CDK6_CRBN_dock*pdb \
#                                             -l test/YKL_06_102_conformer*pdb \
#                                             -da test/decoy_atom_list.txt \
#                                             -la test/linker_atom_list.txt \
#                                             -c 0.4 \
#                                             -r rmsd.txt \
#                                             -t default \
#                                             -wd test/decoy_atom_delete.txt \
#                                             -ld test/linker_atom_delete.txt


# # Test rmsd output
# echo "Test #2"
# python ternary_model_prediction_refactor.py -d test/CDK6_CRBN_dock*pdb \
#                                             -l test/YKL_06_102_conformer*pdb \
#                                             -da test/decoy_atom_list.txt \
#                                             -la test/linker_atom_list.txt \
#                                             -c 0.4 \
#                                             -r rmsd_scocres.txt \
#                                             -t default \
#                                             -wd test/decoy_atom_delete.txt \
#                                             -ld test/linker_atom_delete.txt

# Test without "delete atoms" and "specify"
# echo "Test #3"
# python ternary_model_prediction_refactor.py -d test/CDK6_CRBN_dock*pdb \
#                                             -l test/YKL_06_102_conformer*pdb \
#                                             -da test/decoy_atom_list.txt \
#                                             -la test/linker_atom_list.txt \
#                                             -c 0.4 \
#                                             -r rmsd_score.txt \
#                                             -t specify \


# # Test default naming
# echo "Test #4"
# python ternary_model_prediction_refactor.py -d test/CDK6_CRBN_dock*pdb \
#                                             -l test/YKL_06_102_conformer*pdb \
#                                             -da test/decoy_atom_list.txt \
#                                             -la test/linker_atom_list.txt \
#                                             -c 10.0 \
#                                             -r rmsd_score.txt \
#                                             -t default \


# Test default naming
echo "Test #5"
python ternary_model_prediction_refactor.py -d test/CDK6_CRBN_dock*pdb \
                                            -l test/YKL_06_102_conformer*pdb \
                                            -da test/decoy_atom_list.txt \
                                            -la test/linker_atom_list.txt \
                                            -c 10.0 \
                                            -r rmsd_score.txt \
                                            -t default \

