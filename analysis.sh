#!/bin/bash
head -n 1 Run-configuration-0.txt
#tail -n 1 Run-configuration-*.txt
for file in $(ls|grep ^Run)
do
tail -n 1 $file
done
