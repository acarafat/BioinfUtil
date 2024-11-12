echo "change_gbk_origin"

# Change origin position of pTi_GV3101.gbk to repC (position 23305..24624)
bioinfutils change_gbk_origin --origin 23304 --input test/dataset/pTi_GV3101.gbk --output test/dataset/pTi_GV3101_repC.gbk


echo "filter_fasta"