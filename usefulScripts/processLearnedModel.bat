#!/bin/bash -l
#insert use commands here
use .ruby-2.2.2
use R-3.4
use .cairo-1.14.2
unuse Python-2.7
use Anaconda3
source activate TensorflowEnv
#source activate TensorflowEnvRH7
#  $1 is the config file, of the form:
#MODEL="modelPrefix" # comments
#indcluding:
#MODEL="model output prefix"
#MODELPARAMS="shared model params"
#TFALIASFILE="File containing ID\tALIAS per line for making logo HTML files"
#RUNTIMEPARAMS=" -b 256 -t 2 "
#TESTDATA="OHC test data"
#TESTPARAMS="test-specific model parameters"
#TESTOUT="a name for the output test predictions"
set -e
LogoGenerator="/ahg/regevdata/users/cgdeboer/Programs1/REDUCE_Suite/bin/LogoGenerator"

source $1


mkdir -p $MODEL/Model_Eval
cat $1 > $MODEL/Model_Eval/conf.sh 
cat $MODEL.params | awk '($1!=-1){print}' > $MODEL/Model_Eval/Params.txt

if [ -z ${TFALIASFILE+x} ]
then
	echo TFALIASFILE not provided
	cat $MODEL/Model_Eval/Params.txt > $MODEL/Model_Eval/Params.named.txt
else
	if [ -e $TFALIASFILE ] 
	then
		echo "Making TF Alias file: "$MODEL/Model_Eval/Params.named.txt
		echo "MotifID	Alias" > $MODEL/Model_Eval/Aliases.txt
		cat $TFALIASFILE >> $MODEL/Model_Eval/Aliases.txt
		paste $MODEL/Model_Eval/Aliases.txt $MODEL/Model_Eval/Params.txt > $MODEL/Model_Eval/Params.named.txt
	else
		echo TFALIASFILE: $TFALIASFILE does not exist
	fi
fi



echo Processing learned parameters for $MODEL - SKIPPING!!!L
if [ -e REMOVETHISWHENIWANTITTOWORK$MODEL.pkdms ] 
then
	echo Processing motif logos...
	mkdir -p $MODEL/Model_Eval/Motifs
	optimizedPKdMsToMatrices.py -i $MODEL.pkdms -o $MODEL/Model_Eval/Motifs -in $MODEL/Model_Eval/Params.named.txt -d
	cat $MODEL/Model_Eval/Params.named.txt |awk '(NR>1){print $1}' | xargs -I {} pwmToREDUCE.py -i $MODEL/Model_Eval/Motifs/{}.pkdm -o $MODEL/Model_Eval/Motifs/{}.xml -n
	cat $MODEL/Model_Eval/Params.named.txt |awk '(NR>1){print $1}' | xargs -I {} $LogoGenerator -file=$MODEL/Model_Eval/Motifs/{}.xml -style=raw -output=$MODEL/Model_Eval/Motifs -height=5
	#PFMs
	cat $MODEL/Model_Eval/Params.named.txt |awk '(NR>1){print $1}' | xargs -I {} pwmToREDUCE.py -i $MODEL/Model_Eval/Motifs/{}.pkdm -o $MODEL/Model_Eval/Motifs/{}.pfm.xml -f
	cat $MODEL/Model_Eval/Params.named.txt |awk '(NR>1){print $1}' | xargs -I {} $LogoGenerator -file=$MODEL/Model_Eval/Motifs/{}.pfm.xml -style=bits_info -output=$MODEL/Model_Eval/Motifs  -height=5
	makeMotifComparisonTable.py -i $MODEL/Model_Eval/Params.named.txt -ip Motifs/  -is .png -o $MODEL/Model_Eval/pkdm_motif_table.html
	makeMotifComparisonTable.py -i $MODEL/Model_Eval/Params.named.txt -ip Motifs/  -is .pfm.png -o $MODEL/Model_Eval/pfm_motif_table.html
else
	echo No PKdM file for $MODEL - skipping logo processing.
fi



if [ -e $MODEL.params ] 
then
	# use $MODEL/Model_Eval/Params.named.txt
	echo Doing something with learned activities
else
	echo No ___ file for $MODEL.
fi


for ((i=0;i<${#TESTDATA[@]};++i)); do
	curTESTDATA="${TESTDATA[i]}"
	curTESTOUT="${TESTOUT[i]}"
	if [ -e $curTESTDATA ] #makes sure the stuff after "then" is only run if the job was not completed last time.
	then
		echo Evaluating $MODEL on $curTESTDATA
		mkdir -p $MODEL/Model_Eval/Predictions/
		if [ ! -e $MODEL/Model_Eval/Predictions/$curTESTOUT.done ] 
		then
			echo predictThermodynamicEnhancosomeModel.py -i  $curTESTDATA -o $MODEL/Model_Eval/Predictions/$curTESTOUT  -v -v -v  -M $MODEL.Model $MODELPARAMS $RUNTIMEPARAMS $TESTPARAMS
			predictThermodynamicEnhancosomeModel.py -i  $curTESTDATA -o $MODEL/Model_Eval/Predictions/$curTESTOUT  -v -v -v  -M $MODEL.Model $MODELPARAMS $RUNTIMEPARAMS $TESTPARAMS
			touch $MODEL/Model_Eval/Predictions/$curTESTOUT.done
		else
			echo already predicted on test data
		fi
		correlationPlot.R $MODEL/Model_Eval/Predictions/$curTESTOUT $MODEL/Model_Eval/Predictions/predicted_vs_actual.$curTESTOUT 	predicted actual
	else
		echo No test data provided: $curTESTDATA
	fi
done
