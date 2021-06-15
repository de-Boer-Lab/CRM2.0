# Regulatory Complexity calculation and scripts

Unfortunately, we have not yet formalized this into an easy-to-run script, and so it appears to be somewhat hacky. Luckily, it's a fairly simple calculation, so is relatively straightforward. 

## What you need
* A model (already trained)
* A list of TFs as they have been used in the model (in a file, one per line)
* A file containing one-hot encoded DNA sequences of the appropriate length -  what you are calculating regulatory complexity on.

## Overview
Here is a breif overview of the steps needed to calculate regulatory complexity:

1. Calculate predicted expression for each sequence using the model
2. For each TF, Calculate predicted expression for each sequence using the model with that TF's concentration set to 0
3. Calculate the regulatory interaction strengths for each TF (change in expression relative to the complete unadultered model) and for each sequence
4. Calculate the Gini coeficient of the regulatory interaction strengths for each TF
5. Calculate regulatory complexity: 1-Gini

## Detailed instructions
### Calculate predicted expression for each sequence using the model
Here is an example:
```bash
predictThermodynamicEnhancosomeModel.py -i  all.EDV.OHC.seq.gz  -o all.EDV.SCUra.preds.gz  -v -v -v  -M ../EBound2_progressive_learning_pTpA_SC-Ura.ACPMBOL.Model  -sl 110  -b 256 -t 2
```
The meanings of the parameters are explained in further detail below.

### Predict expression, setting each [TF] to 0 (TF dropout)
The following submits a job to the cluster (ours is based on SGE), which runs the prediction with TF dropout for each TF
```bash
qsub  -b y -cwd  -o dropoutTFsOneByOne_SCUra.out -l h_vmem=50g,h_rt=99:00:00 -e dropoutTFsOneByOne_SCUra.out -q regev -pe smp 2 -N '..dropoutTFsOneByOne.bat .home.unix.cgde' -t 1-245 './dropoutTFsOneByOne.bat /home/unix/cgdeboer/CIS-BP/YeTFaSCo/allTF_PKdMFiles_polyA_and_FZF1.txt 20181110_TF_Dropout/all.SCUra.out EBound2_progressive_learning_pTpA_SC-Ura.ACPMBOL all.EDV.OHC.seq.gz'
```

Within `dropoutTFsOneByOne.sh`, this is the most important line and you just need to run this once per TF:
```bash
predictThermodynamicEnhancosomeModel.py -i /home/unix/cgdeboer/SingleCellPGM/GPRA/20180307_pTpA_80bp_Twist/justSeqs_native_HQ.sequencedRegion110.OHC.gz -res $3.ckpt -o $2.$3.native_HQ.drop.$id.out.gz  -v -v -v  -t 2 -b 256  -sl 110 -M $3.Model -eb -dm /home/unix/cgdeboer/CIS-BP/YeTFaSCo/allTF_PKdMFiles_polyA_and_FZF1.txt -dotf $id $4 $5 $6 $7 $8 $9
# Here are the parameters broken down:
# input sequences: -i <input_file_OHC_sequences>
# Model restore point: -res $3.ckpt 
# Output file: -o $2.$3.native_HQ.drop.$id.out.gz  
# Verbose output: -v -v -v  
# Threads: -t 2 
# Minibatch size: -b 256  
# Length of the sequences: -sl 110 
# Model save file: -M $3.Model 
# Model config; expected binding: -eb 
# Motif file: -dm /home/unix/cgdeboer/CIS-BP/YeTFaSCo/allTF_PKdMFiles_polyA_and_FZF1.txt 
# WHICH TF TO DROPOUT: -dotf <TF_i> 
# Other parameters as given to dropoutTFsOneByOne.sh: $4 $5 $6 $7 $8 $9
```
Here, `<TF_i>` is a value from 0..(#TFs -1). Once complete, these can be concatenated together like so:
```bash
paste 20181110_TF_Dropout/all.SCUra.out.EBound2_progressive_learning_pTpA_SC-Ura.ACPMBOL.*.out | awk '{ for (i=3;i<=NF;i+=2) $i="" } 1' > all.EDV.Gini.dropout.predicitons.SCUra.txt
```

### Calculate regulatory complexity
Here is an example of code used to calculate regulatory complexity in R. 
```R
library(ineq)

#read in the predicted expression levels when each TF is (separately) dropped out of the model
neDropoutPredictions = read.table(sprintf("%s/../20190118_manualSeqDesign/Gini.dropout.predicitons.txt",inDir),header=T,stringsAsFactors=F, na.strings = c("NA","NaN"));
neDropoutPredictions$actual=NULL; #wasn't actually used for anything
names(neDropoutPredictions)=(1:ncol(neDropoutPredictions))-1; #rename columns to the TF indeces

#read in the predicted expression for the complete model
neStandardPredictions = read.table(sprintf("%s/../20190118_manualSeqDesign/pTpA_Glu.ACPMB.pos.Gini.preds.gz",inDir),header=T,sep="\t",stringsAsFactors=F, na.strings = c("NA","NaN"));

#evolvedSeqs is a data.frame containing additional information about the sequences we were predicting on. The order of the sequences in evolvedSeqs, neStandardPredictions, and neDropoutPredictions is the same.
neDropoutPredictions = cbind(cbind(evolvedSeqs, neStandardPredictions["predicted"]), neDropoutPredictions)
neDropoutPredictionsMelted = melt(neDropoutPredictions, id.vars=c("pid","N80seq","giniType","sourcePID","editDistance","predicted"))
names(neDropoutPredictionsMelted)[ncol(neDropoutPredictionsMelted) - (1:0)] = c("i","doPred") #doPred= drop out predicted expression
neDropoutPredictionsMelted$i = as.numeric(as.character(neDropoutPredictionsMelted$i)) # these were column names as strings; convert to int

#Calculate the regulatory interaction strengths for each TF/promoter sequence
neDropoutPredictionsMelted$delOverWT = neDropoutPredictionsMelted$doPred - neDropoutPredictionsMelted$predicted;
#add information about TFs to the regulatory interaction strengths; allTFs is a data.frame with column i the TF index, 
#  and other information about each TF in other columns
neDropoutPredictionsMelted = merge(neDropoutPredictionsMelted, allTFs, by="i")

# Calculate the Gini coeficient for regulatory interaction strength inequality for each sequence
# Here, I have excluded TFs that were not actually contributing to the model (unusedTFs); these can be identified as those which always have regulatory interaction strengths near 0. Even if these non-contributing TFs are included here, the relative values of Gini and regulatory complexity should be unaffected (since they would contribute nothing for all promoters equally).
neDropoutGinis = cast(neDropoutPredictionsMelted[!(neDropoutPredictionsMelted$TF %in% unusedTFs),], sourcePID + editDistance+giniType +predicted~ ., value="delOverWT", fun.aggregate = function(x){ineq(abs(x),type="Gini")})
names(neDropoutGinis)[ncol(neDropoutGinis)]="gini"

#calculate regulatory complexity for each sequence
neDropoutGinis$regulatoryComplexity = 1- neDropoutGinis$gini;
```

