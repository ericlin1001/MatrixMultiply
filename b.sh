#!/bin/bash
dir=../Run/$RANDOM
logFile=runLog.log
rm a.out
mpic++ -DUSE_MPI main.cpp
echo Run:$dir
mkdir $dir
cp -r . $dir/
echo sleep 1
sleep 1
echo Start Run...........................
cd $dir
echo "mpirun -f machinefile -n $1 ./a.out|tee $logFile" > RunCmd.txt
mpirun -f machinefile -n $1 ./a.out|tee $logFile
echo "Results are saved in $dir"

