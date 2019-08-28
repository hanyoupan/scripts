#!/bin/bash                                                                                                                                                                       

cd /nics/c/home/hanyou/scripts/rat_genome_03

for i in 34 39
do
bash 2019-06-13_2_supernova_run_mkoutput_qsub_generator.sh $i
done
