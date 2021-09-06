#!/bin/bash
#$1 = training data
#$2 = name prefix
#$3 = seq length
#$4 ... = other params
set -e #make it so that the script exits completely if one command fails

export curFile=$2.A.done
if [ ! -e $curFile ] 
then
	echo makeThermodynamicEnhancosomeModel.py -i $1 -o ./$2.A -sl $3 -b 1024  -r 1 -l1 0.00001 -l2 0.000001 -t 2 -nm 245 -ml 25 -dm /home/unix/cgdeboer/CIS-BP/YeTFaSCo/allTF_PKdMFiles_polyA_and_FZF1.txt -se 5000 -lr 0.04 -ntc -ntm  -ic /home/unix/cgdeboer/CIS-BP/YeTFaSCo/allTF_minKds_polyA_and_FZF1_justKds.txt $4 $5 $6 $7 $8 $9
	makeThermodynamicEnhancosomeModel.py -i $1 -o ./$2.A -sl $3 -b 1024  -r 1 -l1 0.00001 -l2 0.000001 -t 2 -nm 245 -ml 25 -dm /home/unix/cgdeboer/CIS-BP/YeTFaSCo/allTF_PKdMFiles_polyA_and_FZF1.txt -se 5000 -lr 0.04 -ntc -ntm  -ic /home/unix/cgdeboer/CIS-BP/YeTFaSCo/allTF_minKds_polyA_and_FZF1_justKds.txt $4 $5 $6 $7 $8 $9
	touch $curFile
fi

export curFile=$2.AP.done
if [ ! -e $curFile ] 
then
	echo makeThermodynamicEnhancosomeModel.py -i $1 -o ./$2.AP -sl $3 -b 1024  -r 1 -l1 0.00001 -l2 0.000001 -t 2 -nm 245 -ml 25 -dm /home/unix/cgdeboer/CIS-BP/YeTFaSCo/allTF_PKdMFiles_polyA_and_FZF1.txt -se 5000 -lr 0.04 -ntc -ntm -po -ic /home/unix/cgdeboer/CIS-BP/YeTFaSCo/allTF_minKds_polyA_and_FZF1_justKds.txt $4 $5 $6 $7 $8 $9
	makeThermodynamicEnhancosomeModel.py -i $1 -o ./$2.AP -sl $3 -b 1024  -r 1 -l1 0.00001 -l2 0.000001 -t 2 -nm 245 -ml 25 -dm /home/unix/cgdeboer/CIS-BP/YeTFaSCo/allTF_PKdMFiles_polyA_and_FZF1.txt -se 5000 -lr 0.04 -ntc -ntm -po -ic /home/unix/cgdeboer/CIS-BP/YeTFaSCo/allTF_minKds_polyA_and_FZF1_justKds.txt $4 $5 $6 $7 $8 $9
	touch $curFile
fi

export curFile=$2.ACP.done
if [ ! -e $curFile ] 
then
	echo makeThermodynamicEnhancosomeModel.py -i $1 -o ./$2.ACP -sl $3 -b 1024  -r 1 -l1 0.00001 -l2 0.000001 -t 2 -nm 245 -ml 25 -dm /home/unix/cgdeboer/CIS-BP/YeTFaSCo/allTF_PKdMFiles_polyA_and_FZF1.txt -se 5000 -ntm -po -res $2.AP.ckpt -ic /home/unix/cgdeboer/CIS-BP/YeTFaSCo/allTF_minKds_polyA_and_FZF1_justKds.txt  $4 $5 $6 $7 $8 $9
	makeThermodynamicEnhancosomeModel.py -i $1 -o ./$2.ACP -sl $3 -b 1024  -r 1 -l1 0.00001 -l2 0.000001 -t 2 -nm 245 -ml 25 -dm /home/unix/cgdeboer/CIS-BP/YeTFaSCo/allTF_PKdMFiles_polyA_and_FZF1.txt -se 5000 -ntm -po -res $2.AP.ckpt -ic /home/unix/cgdeboer/CIS-BP/YeTFaSCo/allTF_minKds_polyA_and_FZF1_justKds.txt  $4 $5 $6 $7 $8 $9
	touch $curFile
fi

export curFile=$2.ACPM.done
if [ ! -e $curFile ] 
then
	echo makeThermodynamicEnhancosomeModel.py -i $1 -o ./$2.ACPM -sl $3 -b 1024  -r 1 -l1 0.00001 -l2 0.000001 -t 2 -nm 245 -ml 25 -dm /home/unix/cgdeboer/CIS-BP/YeTFaSCo/allTF_PKdMFiles_polyA_and_FZF1.txt -se 5000 -po -lr 0.001 -res $2.ACP.ckpt  $4 $5 $6 $7 $8 $9 ; touch $curFile
	makeThermodynamicEnhancosomeModel.py -i $1 -o ./$2.ACPM -sl $3 -b 1024  -r 1 -l1 0.00001 -l2 0.000001 -t 2 -nm 245 -ml 25 -dm /home/unix/cgdeboer/CIS-BP/YeTFaSCo/allTF_PKdMFiles_polyA_and_FZF1.txt -se 5000 -po -lr 0.001 -res $2.ACP.ckpt  $4 $5 $6 $7 $8 $9 ; touch $curFile
fi

export curFile=$2.ACPMB.done
if [ ! -e $curFile ] && [ -e $2.ACPM.done ]
then
	echo makeThermodynamicEnhancosomeModel.py -i $1 -o ./$2.ACPMB -sl $3 -b 1024  -r 1 -l1 0.00001 -l2 0.000001 -t 2 -nm 245 -ml 25 -dm /home/unix/cgdeboer/CIS-BP/YeTFaSCo/allTF_PKdMFiles_polyA_and_FZF1.txt -se 5000 -po -lr 0.001 -res $2.ACPM.ckpt -bl -ibl pTpA_Glu.binding_limits.from_elbows.txt $4 $5 $6 $7 $8 $9 ; touch $curFile
	makeThermodynamicEnhancosomeModel.py -i $1 -o ./$2.ACPMB -sl $3 -b 1024  -r 1 -l1 0.00001 -l2 0.000001 -t 2 -nm 245 -ml 25 -dm /home/unix/cgdeboer/CIS-BP/YeTFaSCo/allTF_PKdMFiles_polyA_and_FZF1.txt -se 5000 -po -lr 0.001 -res $2.ACPM.ckpt -bl -ibl pTpA_Glu.binding_limits.from_elbows.txt $4 $5 $6 $7 $8 $9 
	touch $curFile
fi

export curFile=$2.ACPMBOL.done
if [ ! -e $curFile ] 
then
	echo makeThermodynamicEnhancosomeModel.py -i $1 -o ./$2.ACPMBOL -sl $3 -b 1024  -r 1 -l1 0.00001 -l2 0.000001 -t 2 -nm 245 -ml 25 -dm /home/unix/cgdeboer/CIS-BP/YeTFaSCo/allTF_PKdMFiles_polyA_and_FZF1.txt -se 5000 -po -stra -posa -res $2.ACPMB.ckpt -lr 0.0001 $4 $5 $6 $7 $8 $9
	makeThermodynamicEnhancosomeModel.py -i $1 -o ./$2.ACPMBOL -sl $3 -b 1024  -r 1 -l1 0.00001 -l2 0.000001 -t 2 -nm 245 -ml 25 -dm /home/unix/cgdeboer/CIS-BP/YeTFaSCo/allTF_PKdMFiles_polyA_and_FZF1.txt -se 5000 -po -stra -posa -res $2.ACPMB.ckpt -lr 0.0001 $4 $5 $6 $7 $8 $9
	touch $curFile
fi

exit

export curFile=$2.ACPO.done
if [ ! -e $curFile ] 
then
	echo makeThermodynamicEnhancosomeModel.py -i $1 -o ./$2.ACPO -sl $3 -b 1024  -r 1 -l1 0.00001 -l2 0.000001 -t 2 -nm 245 -ml 25 -dm /home/unix/cgdeboer/CIS-BP/YeTFaSCo/allTF_PKdMFiles_polyA_and_FZF1.txt -se 5000 -ntm -po -stra -res $2.ACP.ckpt  $4 $5 $6 $7 $8 $9
	makeThermodynamicEnhancosomeModel.py -i $1 -o ./$2.ACPO -sl $3 -b 1024  -r 1 -l1 0.00001 -l2 0.000001 -t 2 -nm 245 -ml 25 -dm /home/unix/cgdeboer/CIS-BP/YeTFaSCo/allTF_PKdMFiles_polyA_and_FZF1.txt -se 5000 -ntm -po -stra -res $2.ACP.ckpt  $4 $5 $6 $7 $8 $9
	touch $curFile
fi

export curFile=$2.ACPOL.done
if [ ! -e $curFile ] 
then
	echo makeThermodynamicEnhancosomeModel.py -i $1 -o ./$2.ACPOL -sl $3 -b 1024  -r 1 -l1 0.00001 -l2 0.000001 -t 2 -nm 245 -ml 25 -dm /home/unix/cgdeboer/CIS-BP/YeTFaSCo/allTF_PKdMFiles_polyA_and_FZF1.txt -se 5000 -ntm -po -stra -posa -res $2.ACPO.ckpt -lr 0.0001 $4 $5 $6 $7 $8 $9
	makeThermodynamicEnhancosomeModel.py -i $1 -o ./$2.ACPOL -sl $3 -b 1024  -r 1 -l1 0.00001 -l2 0.000001 -t 2 -nm 245 -ml 25 -dm /home/unix/cgdeboer/CIS-BP/YeTFaSCo/allTF_PKdMFiles_polyA_and_FZF1.txt -se 5000 -ntm -po -stra -posa -res $2.ACPO.ckpt -lr 0.0001 $4 $5 $6 $7 $8 $9
	touch $curFile
fi
