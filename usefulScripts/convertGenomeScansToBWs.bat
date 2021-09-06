#!/bin/bash -l
#insert use commands here
#  $1 is the task file

export id=`awk "NR==$SGE_TASK_ID" $1` # input the $SGE_TASK_IDth line of this file into $id
#uncomment and modify the next two lines if your task file has multiple fields.
#splitID=($id) #split id on whitespace into array that can be accessed like ${splitID[0]}
#export id=${splitID[0]} # if there was whitesapce delim fields in the task file, make the new task ID the first field (0)
echo  $SGE_TASK_ID : $id #print the current task to the primary output log
export logDir="20170307_genome_tracks" #fill this in with where you would like the output logs to go - this directory must actually exist.
export logPre=$logDir/$id #a prefix to use for the job output log ($logPre.olog) and the file denoting that the job completed ($logPre.done)
#redirect stdout and stderr to file
exec 1>>$logPre.olog # redirect standard out to an output log file based off of this task ID
exec 2>&1 #redirect stderr to stdout
set -e #make it so that the script exits completely if one command fails


if [ ! -e $logPre.done ] #makes sure the stuff after "then" is only run if the job was not completed last time.
then
	echo Job $JOB_ID:$SGE_TASK_ID started on $HOST: `date` #print to the output file the current job id and task, and host
	#PLACE COMMANDS HERE
	paste <(gzip -dc   temp_20161224_FASTAs/wholeGenomeFasta.110.names.gz | sed 's/-/:/g' | awk 'BEGIN{FS=":";OFS="\t"};{print $1, ($3+$2)/2, ($3+$2)/2+1}')  <(cat 20161129_testing_all_models/predicted.$id.genome.txt | awk 'BEGIN{OFS="\t";FS="\t"};{if (NR>1) {print $3}}')  > 20170307_genome_tracks/$id.genome.pOpen.bedgraph
	bedGraphToBigWig 20170307_genome_tracks/$id.genome.pOpen.bedgraph ~/genomes/sc/20110203_R64/chrom.sizes 20170307_genome_tracks/$id.genome.pOpen.bw
	paste <(gzip -dc   temp_20161224_FASTAs/wholeGenomeFasta.110.names.gz | sed 's/-/:/g' | awk 'BEGIN{FS=":";OFS="\t"};{print $1, ($3+$2)/2, ($3+$2)/2+1}')  <(cat 20161129_testing_all_models/predicted.$id.genome.txt | awk 'BEGIN{OFS="\t";FS="\t"};{if (NR>1) {print $2}}')  > 20170307_genome_tracks/$id.genome.EL.bedgraph
	bedGraphToBigWig 20170307_genome_tracks/$id.genome.EL.bedgraph ~/genomes/sc/20110203_R64/chrom.sizes 20170307_genome_tracks/$id.genome.EL.bw
	rm 20170307_genome_tracks/$id.genome.EL.bedgraph 20170307_genome_tracks/$id.genome.pOpen.bedgraph
	touch $logPre.done #always the last command - create the file that indicates this task completed successfully
	echo Job finished: `date`
else
	echo Job already done. To redo, rm $logPre.done
fi
qstat -j $JOB_ID | grep "^usage *$SGE_TASK_ID:" #display job usage stats
