#!/bin/bash

number=$1

id=$number
time=140:00:00
ppn=24
node_type="partition=monster"

cd /nics/c/home/hanyou/scripts/rat_genome_03

file_name="`date +"%Y-%m-%d"`_x_supernova_run_mkoutput_general_ut_$id.qsub"
bash /nics/c/home/hanyou/scripts/_general/pbs_ut.sh $file_name $time $ppn $node_type > $file_name

echo "source /nics/c/home/hanyou/.bashrc
cd

mkdir -p /lustre/haven/proj/UTHSC0013/UM/projects/rat_genome_03/supernova/$id
cd /lustre/haven/proj/UTHSC0013/UM/projects/rat_genome_03/supernova/$id

supernova run --id=$id \
--fastqs=/lustre/haven/proj/UTHSC0013/UM/projects/rat_genome_03/sequences/$id/ \
--maxreads=1640000000

# maxread = 2.87e9 * 56 / 151 = 1.64e9

mkdir -p fasta
cd fasta

supernova mkoutput \
--style=pseudohap \
--asmdir=/lustre/haven/proj/UTHSC0013/UM/projects/rat_genome_03/supernova/${id}/${id}/outs/assembly \
--outprefix=${id}_pseudohap \
--minsize=10000 \
--headers=short


" >> $file_name

qsub $file_name

