#!/bin/bash

cd /nics/c/home/hanyou/scripts/rat_genome_03

for i in {23..30}
do
bash 2019-06-04_4_long_ranger_qsub_generator.sh $i
done
