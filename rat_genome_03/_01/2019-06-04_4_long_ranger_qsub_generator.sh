#!/bin/bash

number=$1

id="5739-JL-00"$number
time=140:00:00
ppn=40
kind="partition=general,feature=skylake"

cd /nics/c/home/hanyou/scripts/rat_genome_03

file_name="`date +"%Y-%m-%d"`_x_long_ranger_wgs_$id.qsub"
bash /nics/c/home/hanyou/scripts/_pbs.sh $file_name $time $ppn $kind > $file_name

echo "source /nics/c/home/hanyou/.bashrc
cd

mkdir -p /nics/c/home/hanyou/UM/data/long_ranger/rat_genome_03/$id
cd /nics/c/home/hanyou/UM/data/long_ranger/rat_genome_03/$id

longranger wgs --id=$id \
--reference=/nics/c/home/hanyou/UM/resources/references/refdata-rn6.fa.concatenate \
--fastqs=/lustre/haven/proj/UTHSC0013/UM/data/sequences/rat_genome_03/$id/ \
--sample=${id}_1,${id}_2,${id}_3,${id}_4 \
--vcmode=gatk:/lustre/haven/proj/UTHSC0013/UM/resources/apps/gatk_4.0.3.0/gatk-4.0.3.0/gatk-package-4.0.3.0-local.jar
" >> $file_name

qsub $file_name
