#!/bin/bash

number=$1

id=$number
time=140:00:00
ppn=40
kind="partition=general,feature=skylake"

cd /nics/c/home/hanyou/scripts/rat_genome_03

file_name="`date +"%Y-%m-%d"`_x_long_ranger_wgs_$id.qsub"
bash /nics/c/home/hanyou/scripts/_general/pbs_ut.sh $file_name $time $ppn $kind > $file_name

echo "source /nics/c/home/hanyou/.bashrc
cd

mkdir -p /lustre/haven/proj/UTHSC0013/UM/projects/rat_genome_03/long_ranger/$id
cd /lustre/haven/proj/UTHSC0013/UM/projects/rat_genome_03/long_ranger/$id

longranger wgs --id=$id \
--reference=/lustre/haven/proj/UTHSC0013/UM/projects/rat_genome_03/references/refdata-rn6.fa.concatenate \
--fastqs=/lustre/haven/proj/UTHSC0013/UM/projects/rat_genome_03/sequences/$id/ \
--vcmode=gatk:/lustre/haven/proj/UTHSC0013/UM/apps/gatk_4.0.3.0/gatk-4.0.3.0/gatk-package-4.0.3.0-local.jar
" >> $file_name

qsub $file_name

