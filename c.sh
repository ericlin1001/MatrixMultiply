#!/bin/bash
rm a.out
sleep 1
mpic++ -DUSE_MPI main.cpp
echo sleep 1
sleep 1
echo Start Run...........................
mpirun -f machinefile -n $1 ./a.out
mv Run* output/
mv Data* output/
echo End Run...........................

