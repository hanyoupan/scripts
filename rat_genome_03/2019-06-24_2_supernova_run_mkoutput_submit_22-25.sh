#!/bin/bash                                                                                                                                                                       

cd /nics/c/home/hanyou/scripts/rat_genome_03

start=22
end=25

for i in $(seq $start $end)
do
bash 2019-06-13_2_supernova_run_mkoutput_qsub_generator.sh $i
done
