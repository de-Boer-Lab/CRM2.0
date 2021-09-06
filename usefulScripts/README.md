Here are some useful scripts. Dependencies are listed within the scripts. Likely many will not run out of the box. As we try to use them and they break, we will add needed documentation. Mostly, these are just for added convenience and are not necessary to do anything needed.

Many of these commands take the form of a submitted job, which includes the command to be run (e.g. if it starts with `submitlog`)

They are ordered roughly from most useful to least. 

## Make a new model
```bash
use UGES; submitlog -g 1 -m 40 -n 2 -o makeProgressiveBiochemModels_OHC_Scratch_pTpA_SC-Ura_EB2.olog  "./makeProgressiveBiochemicalModels.bat 20161024_average_promoter_ELs_per_seq_N80_SC-Ura_Y8203_ALL.shuffled_OHC_train.txt.gz EBound2_progressive_learning_pTpA_SC-Ura 110 -eb -v -v -v"
```

## Testing a model and processing learned parameters
```bash
use UGES; submitlog -m 40 -n 2 -g 1 -o process.pTpA_SC-Ura.ACPMBOL.olog "./processLearnedModel.bat EBound2_progressive_learning_pTpA_SC-Ura.ACPMBOL.conf.sh"
```


## Regulatory complexity
```bash
use UGES; submitlog -m 40 -n 2 -g 0 -o predictOn1002YeastGenomeOrthopros.SCUra.olog "predictThermodynamicEnhancosomeModel.py -i  20190719_1002_yeast_genomes/orthoNative80.2.110.OHC.seq.gz  -o 20190719_1002_yeast_genomes/orthoNative80.pTpA_SCUra.preds.gz  -v -v -v  -M EBound2_progressive_learning_pTpA_SC-Ura.ACPMBOL.Model  -sl 110  -b 256 -t 2 "
submitlog -n 2 -m 50 -g 0 -t /home/unix/cgdeboer/CIS-BP/YeTFaSCo/allTF_PKdMFiles_polyA_and_FZF1.txt  -o calcGiniByTFDropout_P1.orthopros.SCUra.olog ./calcGiniByTFDropout_P1.bat /home/unix/cgdeboer/CIS-BP/YeTFaSCo/allTF_PKdMFiles_polyA_and_FZF1.txt 20190719_1002_yeast_genomes/20190719_TF_Dropout/all.SCUra.out  EBound2_progressive_learning_pTpA_SC-Ura.ACPMBOL  20190719_1002_yeast_genomes/orthoNative80.2.110.OHC.seq.gz
paste 20190719_1002_yeast_genomes/20190719_TF_Dropout//all.SCUra.out.EBound2_progressive_learning_pTpA_SC-Ura.ACPMBOL.*.out | awk 'BEGIN{OFS="\t"}{ for (i=3;i<=NF;i+=2) $i="" } 1' | sed 's/\t\t/\t/g' | gzip  > 20190719_1002_yeast_genomes/orthopros.SCUra.Gini.dropout.predictions.txt.gz
# at this point, a file has been created (20190719_1002_yeast_genomes/orthopros.SCUra.Gini.dropout.predictions.txt.gz) that contains the regulatory interactions with which to calculate the regulatory complexity
```


## Make bigwig files of genomic scans
```bash
#scan whole genome with models
#make whole genome fasta file
submitlog -m 1 -o makeWholeGenomeFasta.sc.olog "tileGenomeFasta.py -i ~/genomes/sc/20110203_R64/chrX.map -s 110 -p 1 -t -o temp_20161224_FASTAs/wholeGenomeFasta.110.seqs"
cat temp_20161224_FASTAs/wholeGenomeFasta.110.seqs | awk 'BEGIN{OFS="\t";FS="\t"};{print $1}' | gzip -c >temp_20161224_FASTAs/wholeGenomeFasta.110.names.gz
cat temp_20161224_FASTAs/wholeGenomeFasta.110.seqs | awk 'BEGIN{OFS="\t";FS="\t"};{print $2,0}' | gzip -c >temp_20161224_FASTAs/wholeGenomeFasta.110.seqs.gz
submitlog -m 1 -o make_OHC_wholeGenome.olog "seqsToOHC.py -i  temp_20161224_FASTAs/wholeGenomeFasta.110.seqs.gz  -m 110 -o temp_20161224_FASTAs/wholeGenomeFasta.110.OHC.gz"
#rerun scanTestDataWithBiochemicalEnhancosomeModels.bat
ls 20161129_testing_all_models/predicted.*genome* -1 | sed 's/20161129_testing_all_models\/predicted.//g' | sed 's/.genome.txt//g' > all_EBound2_Model_Namess.txt
submitlog -m 5 -o convertGenomeScansToBWs.olog -t all_EBound2_Model_Namess.txt ./convertGenomeScansToBWs.bat all_EBound2_Model_Namess.txt
```

## make PNGs of motifs and HTML of PNGs
```bash
submitlog -m 5 -o motifsToPNG.pTpA_Glu.ACPM.olog ./convertLearnedPWMsToPNGs.bat ../../../20160525_NextSeq_pTpA_and_Abf1TATA/analysis/tensorflow/EBound2_progressive_learning_pTpA_Glu.ACPM.pkdms LEARNED_LOGOS/EBound2_progressive_learning/pTpA_Glu.ACP
#compile all such motifs into an HTML
makeMotifComparisonTable.py -i ~/CIS-BP/YeTFaSCo/allTF_aliases_polyA_and_FZF1.txt -ip Original/,EBound2_progressive_learning/pTpA_Glu.ACPM/,EBound2_progressive_learning/Abf1TATA_Glu.ACPM/,EBound2_progressive_learning/pTpA_Gal.ACPM/,EBound2_progressive_learning/pTpA_Gly.ACPM/ -is .pfm.png,.pfm.png,.pfm.png,.pfm.png,.pfm.png -o LEARNED_LOGOS/pfm_motif_table.html
```

## Activity class model
Experimental models containing different types of TF activities that independently contribute to gene expression
```bash
submitlog -m 40 -n 2 -g 1 -o trainActivityClassesModel.AC4.pTpA_Glu.olog ./trainActivityClassesModel_nAC.bat sisterDir/20160609_average_promoter_ELs_per_seq_pTpA_ALL.shuffled_OHC_train.txt.gz pTpA_Glu 110 4
```

## early version of regulatory complexity
```bash
submitlog -n 2 -m 50 -g 1 -t /home/unix/cgdeboer/CIS-BP/YeTFaSCo/allTF_PKdMFiles_polyA_and_FZF1.txt  -o dropoutTFsOneByOne.out ./dropoutTFsOneByOne.bat /home/unix/cgdeboer/CIS-BP/YeTFaSCo/allTF_PKdMFiles_polyA_and_FZF1.txt 20181110_TF_Dropout/out  EBound2_progressive_learning_pTpA_Glu.ACPMB.pos " -VARIABLE"
```


## convert old-style models to new
This is for upgrading old versions of models. Should not be needed for new models.
```bash
use UGES; submitlog -o convertOldModelsToNewModels.olog -t modelsToMakeFromCkpts.txt -n 2 -g 0 -m 60 ./convertOldModelsToNewModels.bat modelsToMakeFromCkpts.txt
```
