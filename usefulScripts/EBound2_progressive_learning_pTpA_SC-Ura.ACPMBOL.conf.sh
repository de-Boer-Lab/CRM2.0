MODELPARAMS=" -sl 110 "
TFALIASFILE="/home/unix/cgdeboer/CIS-BP/YeTFaSCo/allTF_aliases_polyA_and_FZF1.txt"
RUNTIMEPARAMS=" -b 256 -t 2 "
TESTPARAMS="  "
MODELCOMMAND="makeThermodynamicEnhancosomeModel.py -i 20161024_average_promoter_ELs_per_seq_N80_SC-Ura_Y8203_ALL.shuffled_OHC_train.txt.gz -o ./EBound2_progressive_learning_pTpA_SC-Ura.ACPMBOL -sl 110 -b 1024 -r 1 -l1 0.00001 -l2 0.000001 -t 2 -nm 245 -ml 25 -dm /home/unix/cgdeboer/CIS-BP/YeTFaSCo/allTF_PKdMFiles_polyA_and_FZF1.txt -se 50 -po -stra -posa -res EBound2_progressive_learning_pTpA_SC-Ura.ACPMB.ckpt -lr 0.0001 -eb -v -v -v"
MODEL="EBound2_progressive_learning_pTpA_SC-Ura.ACPMBOL"
TESTDATA=( 20161024_average_promoter_ELs_per_seq_N80_SC-Ura_Y8203_ALL.shuffled_OHC_test.txt.gz /home/unix/cgdeboer/SingleCellPGM/GPRA/20160525_NextSeq_pTpA_and_Abf1TATA/analysis/tensorflow/20160503_average_promoter_ELs_per_seq_atLeast100Counts.OHC.txt.gz  )
TESTOUT=( LQ_test_SC-Ura_pTpA HQ_test_Glu_pTpA )
