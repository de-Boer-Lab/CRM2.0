#!/bin/bash
#$1 = training data
#$2 = id
#$3 = seq len
#$4 = nAC
if [ ! -e EBound3_progressive_learning.$2.AC$4.APC.done ]
then
	makeThermodynamicEnhancosomeModel.py -i $1 -o ./EBound3_progressive_learning.$2.AC$4.APC -sl $3 -b 1024 -r 1 -l1 0.00001 -l2 0.00001 -t 2 -nm 245 -ml 30 -dm /home/unix/cgdeboer/CIS-BP/YeTFaSCo/allTF_PKdMFiles_polyA_and_FZF1.txt -se 10000 -po -eb -v -v -v -raw 1000 -ac $4 -ia 1.0 -ip 0.01 -nta -ntm -ic ~/CIS-BP/YeTFaSCo/allTF_minKds_polyA_and_FZF1_justKds.txt -lr 0.05 -lred 0.98
else
	echo Done APC
fi
if [ ! -e EBound3_progressive_learning.$2.AC$4.APCM.done ]
then
	makeThermodynamicEnhancosomeModel.py -i $1 -o ./EBound3_progressive_learning.$2.AC$4.APCM -sl $3 -b 1024 -r 1 -l1 0.00001 -l2 0.00001 -t 2 -nm 245 -ml 30 -dm /home/unix/cgdeboer/CIS-BP/YeTFaSCo/allTF_PKdMFiles_polyA_and_FZF1.txt -se 10000 -po -eb -v -v -v -raw 1000 -ac $4 -ia 1.0 -nta -ic ~/CIS-BP/YeTFaSCo/allTF_minKds_polyA_and_FZF1_justKds.txt -res EBound3_progressive_learning.$2.AC$4.APC.ckpt -lr 0.01 -lred 0.98 -ca
else
	echo Done APCM
fi
if [ ! -e EBound3_progressive_learning.$2.AC$4.APCM.pos.done ]
then
	makeThermodynamicEnhancosomeModel.py -i $1 -o ./EBound3_progressive_learning.$2.AC$4.APCM.pos -sl $3 -b 1024 -r 1 -l1 0.00001 -l2 0.00001 -t 2 -nm 245 -ml 30 -dm /home/unix/cgdeboer/CIS-BP/YeTFaSCo/allTF_PKdMFiles_polyA_and_FZF1.txt -se 10000 -po -eb -v -v -v -raw 1000 -ac $4 -ia 1.0 -nta -ic ~/CIS-BP/YeTFaSCo/allTF_minKds_polyA_and_FZF1_justKds.txt -res EBound3_progressive_learning.$2.AC$4.APCM.ckpt -lr 0.001 -lred 0.98 -ca -posa -stra
else
	echo Done APCM.pos
fi

