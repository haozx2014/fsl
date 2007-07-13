#!/bin/sh

subjdir=$1
nproc=$2

#$ -S /bin/sh
#$ -V
#$ -N ps_postproc
#$ -m ae

#POSSUMDIR=~ivana/fsldev/${FSLMACHTYPE}/bin
echo $POSSUMDIR

echo Summing all signal from different proccesses into one total signal

${POSSUMDIR}/proc_sum -i ${subjdir}/diff_proc/signal_proc_ \
-o ${subjdir}/signal -n $nproc


echo Converting the signal into the image

${POSSUMDIR}/signal2image -i ${subjdir}/signal \
-o ${subjdir}/image -p ${subjdir}/pulse -a 


echo Removing intermediate files

if [ `imtest ${subjdir}/image_abs` -eq 1 ];then
      rm -rf ${subjdir}/diff_proc
      gzip --best pulse
fi

echo Done
