#!/bin/bash
#$1 = input file of pkdms
#$2 = outDir

mkdir -p $2
optimizedPKdMsToMatrices.py -i $1 -o $2 -in ~/CIS-BP/YeTFaSCo/allTF_aliases_polyA_and_FZF1.txt
cat ~/CIS-BP/YeTFaSCo/allTF_aliases_polyA_and_FZF1.txt |awk '{print $1}' | xargs -I {} pwmToREDUCE.py -i $2/{}.pkdm -o $2/{}.xml -n
cat ~/CIS-BP/YeTFaSCo/allTF_aliases_polyA_and_FZF1.txt |awk '{print $1}' | xargs -I {} LogoGenerator -file=$2/{}.xml -style=raw -output=$2 -height=5
#PFMs
cat ~/CIS-BP/YeTFaSCo/allTF_aliases_polyA_and_FZF1.txt |awk '{print $1}' | xargs -I {} pwmToREDUCE.py -i $2/{}.pkdm -o $2/{}.pfm.xml -f
cat ~/CIS-BP/YeTFaSCo/allTF_aliases_polyA_and_FZF1.txt |awk '{print $1}' | xargs -I {} LogoGenerator -file=$2/{}.pfm.xml -style=bits_info -output=$2  -height=5
