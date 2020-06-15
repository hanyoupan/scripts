#!/bin/bash                                                                                                                                                                       

cd /nics/c/home/hanyou/scripts/rat_genome_03


for i in 25 29 30 32 33 36 37 38 40
do
bash 2019-11-08_1_supernova_run_mkoutput_qsub_generator_20HRLR.sh $i
done
