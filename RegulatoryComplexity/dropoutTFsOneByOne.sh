#!/bin/bash -l
#insert use commands here
#  $1 is the task file
# $2 is the output prefix
# $3 is the model
# $4..$9 are any additional params needed (e.g. -VARIABLE)
export id=`expr $SGE_TASK_ID - 1 `
#uncomment and modify the next two lines if your task file has multiple fields.
#splitID=($id) #split id on whitespace into array that can be accessed like ${splitID[0]}
#export id=${splitID[0]} # if there was whitesapce delim fields in the task file, make the new task ID the first field (0)
#export pfm=${splitID[0]} # if there was whitesapce delim fields in the task file, make the new task ID the first field (0)
echo  $SGE_TASK_ID : $id #print the current task to the primary output log
export logDir="20181110_TF_Dropout" #fill this in with where you would like the output logs to go - this directory must actually exist.
export logPre=$logDir/$id #a prefix to use for the job output log ($logPre.olog) and the file denoting that the job completed ($logPre.done)
#redirect stdout and stderr to file
exec 1>>$logPre.olog # redirect standard out to an output log file based off of this task ID
exec 2>&1 #redirect stderr to stdout
set -e #make it so that the script exits completely if one command fails


if [ ! -e $logPre.done ] #makes sure the stuff after "then" is only run if the job was not completed last time.
then
  echo Job $JOB_ID:$SGE_TASK_ID started on $HOST: `date` #print to the output file the current job id and task, and host
  #PLACE COMMANDS HERE
  predictThermodynamicEnhancosomeModel.py -i /home/unix/cgdeboer/SingleCellPGM/GPRA/20180307_pTpA_80bp_Twist/justSeqs_native_HQ.sequencedRegion110.OHC.gz -res $3.ckpt -o $2.$3.native_HQ.drop.$id.out.gz  -v -v -v  -t 2 -b 256  -sl 110 -M $3.Model -eb -dm /home/unix/cgdeboer/CIS-BP/YeTFaSCo/allTF_PKdMFiles_polyA_and_FZF1.txt -dotf $id $4 $5 $6 $7 $8 $9
  predictThermodynamicEnhancosomeModel.py -i /home/unix/cgdeboer/SingleCellPGM/GPRA/20160525_NextSeq_pTpA_and_Abf1TATA/analysis/tensorflow/20160503_average_promoter_ELs_per_seq_atLeast100Counts.OHC.txt.gz -res $3.ckpt -o $2.$3.random_HQ.drop.$id.out.gz  -v -v -v  -t 2 -b 256  -sl 110 -M $3.Model -eb -dm /home/unix/cgdeboer/CIS-BP/YeTFaSCo/allTF_PKdMFiles_polyA_and_FZF1.txt -dotf $id $4 $5 $6 $7 $8 $9
  touch $logPre.done #always the last command - create the file that indicates this task completed successfully
  echo Job finished: `date`
else
  echo Job already done. To redo, rm $logPre.done
fi
qstat -j $JOB_ID | grep "^usage *$SGE_TASK_ID:" #display job usage stats
