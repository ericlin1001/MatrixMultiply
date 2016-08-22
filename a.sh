#!/bin/bash
dir=../Run/$RANDOM
echo Run:$dir
mkdir $dir
cp -r . $dir/
sleep 1
echo Start Run...........................
mpirun -f machinefile -n 10 $dir/a.out
#mpirun -n 4 $dir/a.out
#rm -rf $dir/*

