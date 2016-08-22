#!/bin/bash
#Edited by EricLin
#Finished:Thu Dec  3 22:31:41 CST 2015
filename=machinefile
echo "checking Machinefile:$filename........."
okFile=${filename}_ok
noFile=${filename}_fail
cat /dev/null>$okFile
cat /dev/null>$noFile
for line in $(cat $filename)
do 
#echo checking machine:$line...;
name=`echo $line|awk -F : '{print $1}'`
ssh cluster@$name "hostname";
if [ "$?" -eq "0" ]; then
echo $line>>$okFile
else
echo $line>>$noFile
fi
done
echo Sucessful lists are saved in file:$okFile
echo Failed lists are saved in file:$noFile
