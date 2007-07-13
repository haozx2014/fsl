#!/bin/sh

subjdir=$1
comm=$2
taskid=`expr $SGE_TASK_ID - 1`

#$ -S /bin/sh
#$ -V -N possum
#$ -m a

#POSSUMDIR=~ivana/fsldev/${FSLMACHTYPE}/bin
echo $POSSUMDIR

echo "${POSSUMDIR}/possum $comm --procid=${taskid} -o ${subjdir}/diff_proc/signal_proc_${taskid}" >  ${subjdir}/diff_proc/script_proc_${taskid}.sh


${POSSUMDIR}/possum $comm --procid=${taskid}  -o ${subjdir}/diff_proc/signal_proc_${taskid}
