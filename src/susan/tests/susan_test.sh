#!/bin/sh
#testing susan
flag=0
#susan 2D gaussian test
susan_smooth  $FSLTESTDIR/common/vol0000 2500 temp 0.01 3 1 1 $FSLTESTDIR/common/vol0090 2400
../susan      $FSLTESTDIR/common/vol0000 2500 0.01 3 1 1 $FSLTESTDIR/common/vol0090 2400 temp1
avwmaths temp -sub temp1 -abs temp2
output1=`avwstats temp2 -m`
echo "Result of susan 3D 3X3x3 1 usan test:" $output1
limit=3  #Good value, test seems sensitive to changes in algorithm/input errors
output2=`echo $output1 $limit | awk '{ print ($1 < $2) ? "passed" : "failed" }' `
echo "$output2"
if [ $output2 != passed ];
then
flag=1
echo "Possible problem with 3D 3x3x3 1 usan test in susan";
else echo "susan 3D 3x3x3 1 usan test within tolerance";
fi
rm temp*

../susan      $FSLTESTDIR/common/vol0000 2500 0.01 3 1 0 temp
../susan      $FSLTESTDIR/common/vol0000 2500 0.01 3 1 1 $FSLTESTDIR/common/vol0000 2500 temp1
avwmaths temp -sub temp1 -abs temp2
output1=`avwstats temp2 -m`
echo "Result of susan 3D 3X3x3 sim usan test:" $output1
limit=3  #Good value, test seems sensitive to changes in algorithm/input errors
output2=`echo $output1 $limit | awk '{ print ($1 < $2) ? "passed" : "failed" }' `
echo "$output2"
if [ $output2 != passed ];
then
flag=1
echo "Possible problem with 3D 3x3x3 sim usan test in susan";
else echo "susan 3D 3x3x3 sim usan test within tolerance";
fi
rm temp*


susan_smooth  $FSLTESTDIR/common/vol0000 2500 temp 0.01 3 1 2 $FSLTESTDIR/common/vol0090 2400 $FSLTESTDIR/common/vol0010 2700
../susan      $FSLTESTDIR/common/vol0000 2500 0.01 3 1 2 $FSLTESTDIR/common/vol0090 2400 $FSLTESTDIR/common/vol0010 2700 temp1
avwmaths temp -sub temp1 -abs temp2
output1=`avwstats temp2 -m`
echo "Result of susan 3D 3X3x3 2 usan test:" $output1
limit=3  #Good value, test seems sensitive to changes in algorithm/input errors
output2=`echo $output1 $limit | awk '{ print ($1 < $2) ? "passed" : "failed" }' `
echo "$output2"
if [ $output2 != passed ];
then
flag=1
echo "Possible problem with 3D 3x3x3 2 usan test in susan";
else echo "susan 3D 3x3x3 2 usan test within tolerance";
fi
rm temp*


susan_smooth  $FSLTESTDIR/common/vol0000 1500 temp 12 2 1 0
../susan      $FSLTESTDIR/common/vol0000 1500 12 2 1 0 temp1
avwmaths temp -sub temp1 -abs temp2
output1=`avwstats temp2 -m`
echo "Result of susan 2D gaussian test:" $output1
limit=3  #Good value, test seems sensitive to changes in algorithm/input errors
output2=`echo $output1 $limit | awk '{ print ($1 < $2) ? "passed" : "failed" }' `
echo "$output2"
if [ $output2 != passed ];
then
flag=1
echo "Possible problem with 2D gaussian test in susan";
else echo "susan 2D gaussian test within tolerance";
fi
rm temp*

susan_smooth  $FSLTESTDIR/common/vol0000 1500 temp 12 3 1 0
../susan      $FSLTESTDIR/common/vol0000 1500 12 3 1 0 temp1
avwmaths temp -sub temp1 -abs temp2
output1=`avwstats temp2 -m`
echo "Result of susan 3D gaussian test:" $output1
limit=4  #Good value, test seems sensitive to changes in algorithm/input errors
output2=`echo $output1 $limit | awk '{ print ($1 < $2) ? "passed" : "failed" }' `
echo "$output2"
if [ $output2 != passed ];
then
flag=1
echo "Possible problem with 3D gaussian test in susan";
else echo "susan 3D gaussian test within tolerance";
fi
rm temp*


susan_smooth  $FSLTESTDIR/common/vol0000 2500 temp 0.01 2 1 0
../susan      $FSLTESTDIR/common/vol0000 2500 0.01 2 1 0 temp1
avwmaths temp -sub temp1 -abs temp2
output1=`avwstats temp2 -m`
echo "Result of susan 2D 3X3 test:" $output1
limit=3  #Good value, test seems sensitive to changes in algorithm/input errors
output2=`echo $output1 $limit | awk '{ print ($1 < $2) ? "passed" : "failed" }' `
echo "$output2"
if [ $output2 != passed ];
then
flag=1
echo "Possible problem with 2D 3x3 test in susan";
else echo "susan 2D 3x3 test within tolerance";
fi
rm temp*

susan_smooth  $FSLTESTDIR/common/vol0000 2500 temp 0.01 3 1 0
../susan      $FSLTESTDIR/common/vol0000 2500 0.01 3 1 0 temp1
avwmaths temp -sub temp1 -abs temp2
output1=`avwstats temp2 -m`
echo "Result of susan 2D 3X3 test:" $output1
limit=3  #Good value, test seems sensitive to changes in algorithm/input errors
output2=`echo $output1 $limit | awk '{ print ($1 < $2) ? "passed" : "failed" }' `
echo "$output2"
if [ $output2 != passed ];
then
flag=1
echo "Possible problem with 2D 3x3 test in susan";
else echo "susan 2D 3x3 test within tolerance";
fi
rm temp*



#../avwmaths++  $FSLTESTDIR/common/filtered_func_data.nii.gz -exp infvol -odt float
#for loop in -nan -nanm 
#do
#avwmaths_32R infvol $loop temp
#../avwmaths++ infvol $loop temp1 -odt float
#avwmaths temp1 -Tmean temp2 
#avwmaths temp -Tmean temp1
#avwmaths temp1 -sub temp2 -abs temp
#output1=`avwstats temp -m`
#echo "Result of $loop volume test:" $output1
#if [ $output1 != 0.000000 ];
#then
#flag=1
#echo "Possible problem with $loop volume test in avwmaths++";
#fi
#rm temp*
#done
#rm infvol*





echo ""
if [ $flag -ne 0 ]
then
echo "non-zero comparison; possible problem with avwmaths++"
else 
echo "All comparisons zero. No problems detected with avwmaths++"
fi