#!/bin/sh
flag=0

for setting in -t -x -y -z
do
../avwsplit++ $FSLTESTDIR/common/filtered_func_data.nii.gz vol $setting
../avwmerge $setting temp vol*
avwmaths temp -Tmean temp2 
avwmaths $FSLTESTDIR/common/filtered_func_data.nii.gz -Tmean temp
avwmaths temp -sub temp2 temp3
output1=`avwstats temp3 -m`
echo "Result of $setting split test:" $output1
if [ $output1 != 0.000000 ];
then
flag=1
echo "Possible problem with $setting test in avwsplit++";
fi
rm temp*
rm vol*
done

echo ""
if [ $flag -ne 0 ]
then
echo "non-zero comparison; possible problem with avwsplit"
else 
echo "All comparisons zero. No problems detected with avwsplit"
fi