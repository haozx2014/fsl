#!/bin/sh

subjdir=$1
nproc=$2

#$ -S /bin/sh
#$ -V
#$ -N ps_postproc
#$ -m ae

run(){
 echo "$1" >> $2/possum.log
 date >> $2/possum.log
 $1 >> $2/possum.log 2>&1
 date >> $2/possum.log
}

echo Summing all signal from different proccesses into one total signal
run "${FSLDIR}/bin/possum_sum -i ${subjdir}/diff_proc/signal_proc_ -o ${subjdir}/signal -n $nproc " ${subjdir}

echo Converting the signal into the image
run "${FSLDIR}/bin/signal2image -i ${subjdir}/signal -o ${subjdir}/image -p ${subjdir}/pulse -a " ${subjdir}

echo Removing intermediate files
if [ `imtest ${subjdir}/image_abs` -eq 1 ];then
      rm -rf ${subjdir}/diff_proc
fi

echo Adding noise
n=`cat ${subjdir}/noise | awk '{print $1 }'`
m=`cat ${subjdir}/noise | awk '{print $2 }'`
Run "echo n $n m $m" $subjdir
fslmaths ${subjdir}/image_abs -Tmean ${subjdir}/image_mean
P98=`fslstats ${subjdir}/image_mean -P 98`
P02=`fslstats ${subjdir}/image_mean -P 2`
tresh=`echo "0.1 * $P98 + 0.9 * $P02 "|bc -l`
fslmaths ${subjdir}/image_mean -thr $tresh ${subjdir}/image_mean
medint=`fslstats ${subjdir}/image_mean -P 50`
dim1=`fslval ${subjdir}/image_abs dim1`
if [ $n == "snr" ]; then
  run "echo Entered the loop 1" $subjdir
  snr=$m
  if [ $snr != 0 ]; then
     sigma=`echo " ${medint} / ( 2 * ${dim1} * ${snr} ) "| bc -l` #I worked this out ages ago.
     echo "sigma ${sigma} snr ${snr} medintensity ${medint}" > ${subjdir}/noise 
  else
     echo "snr  0" > ${subjdir}/noise
  fi
else
  run "echo Entered the loop 2" $subjdir 
  sigma=$m
  if [ $sigma != 0]; then
     snr=`echo " ${medint} / ( 2 * ${dim1} * ${sigma} ) "| bc -l`
     echo "sigma $sigma snr $snr" > ${subjdir}/noise
  fi
fi
run "${FSLDIR}/bin/systemnoise --in=${subjdir}/signal --out=${subjdir}/signalSNR${snr} --sigma=${sigma} " $subjdir
run "${FSLDIR}/bin/signal2image -i ${subjdir}/signalSNR${snr} -o ${subjdir}/image -p ${subjdir}/pulse -a " $subjdir
rm ${subjdir}/signalSNR$snr 
rm ${subjdir}/image_mean.nii.gz

echo Done
