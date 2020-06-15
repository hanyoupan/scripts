#!/bin/bash

cd /nics/c/home/hanyou/scripts/rat_genome_03

for i in {21..40}
do
bash 2020-05-21_1_long_ranger_qsub_generator_nor1.sh "5739-JL-00"$i
done
