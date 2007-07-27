#!/bin/sh

#$ -cwd
#$ -q short.q
#$ -S /bin/sh
#$ -V
#$ -N bpx_preproc
#$ -m as

subjdir=$1

echo Copying files to bedpost directory
cp ${subjdir}/bvecs ${subjdir}/bvals ${subjdir}.bedpostX
${FSLDIR}/bin/imcp ${subjdir}/nodif ${subjdir}/nodif_brain_mask ${subjdir}.bedpostX
${FSLDIR}/bin/avwmaths\
 ${subjdir}/nodif\
 -mas ${subjdir}/nodif_brain_mask\
 ${subjdir}.bedpostX/nodif_brain

${FSLDIR}/bin/avwslice ${subjdir}/data
${FSLDIR}/bin/avwslice ${subjdir}/nodif_brain_mask
echo Done