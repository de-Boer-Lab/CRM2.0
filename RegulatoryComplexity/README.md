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
```{bash}
predictThermodynamicEnhancosomeModel.py -i /home/unix/cgdeboer/SingleCellPGM/GPRA/20180307_pTpA_80bp_Twist/justSeqs_native_HQ.sequencedRegion110.OHC.gz -res $3.ckpt -o $2.$3.native_HQ.drop.$id.out.gz  -v -v -v  -t 2 -b 256  -sl 110 -M $3.Model -eb -dm /home/unix/cgdeboer/CIS-BP/YeTFaSCo/allTF_PKdMFiles_polyA_and_FZF1.txt
```
The meanings of the parameters are explained in further detail below.

### Predict expression, setting each [TF] to 0 (TF dropout)
The following submits a job to the cluster (ours is based on SGE), which runs the prediction with TF dropout for each TF
```{bash}
qsub  -b y -cwd  -o dropoutTFsOneByOne.out -l h_vmem=20g,h_rt=99:00:00 -e dropoutTFsOneByOne.out -q regev -l gpu=1 -pe smp 2 -N '..dropoutTFsOneByOne.sh .home.unix.cgde' -t 1-245 './dropoutTFsOneByOne.sh /home/unix/cgdeboer/CIS-BP/YeTFaSCo/allTF_PKdMFiles_polyA_and_FZF1.txt 20181110_TF_Dropout/out EBound2_progressive_learning_pTpA_Glu.ACPMB.pos  -VARIABLE'
```

Within `dropoutTFsOneByOne.sh`, this is the most important line and you just need to run this once per TF:
```{bash}
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
Here, `<TF_i>` is a value from 0..(#TFs -1).

### Calculate regulatory complexity
In [this file](RegComp.R) I have provided an example of code used to calculate regulatory complexity in R. 


