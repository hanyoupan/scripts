#!/bin/bash

number=$1

id="5739-JL-00"$number
time=140:00:00
ppn=24
node_type="partition=monster"

cd /nics/c/home/hanyou/scripts/rat_genome_03

file_name="`date +"%Y-%m-%d"`_x_supernova_run_mkoutput_$id.qsub"
bash /nics/c/home/hanyou/scripts/_pbs.sh $file_name $time $ppn $node_type > $file_name

echo "source /nics/c/home/hanyou/.bashrc
cd

mkdir -p /lustre/haven/proj/UTHSC0013/UM/data/supernova/rat_genome_03/$id
cd /lustre/haven/proj/UTHSC0013/UM/data/supernova/rat_genome_03/$id

supernova run --id=$id \
--fastqs=/lustre/haven/proj/UTHSC0013/UM/data/sequences/rat_genome_03/${id}/ \
--sample=${id}_1,${id}_2,${id}_3,${id}_4 \
--maxreads=1640000000

# maxread = 2.87e9 * 56 / 151 = 1.64e9

cd /lustre/haven/proj/UTHSC0013/UM/data/supernova/rat_genome_03/$id
mkdir -p fasta
cd fasta

supernova mkoutput \
--style=pseudohap \
--asmdir=/lustre/haven/proj/UTHSC0013/UM/data/supernova/rat_genome_03/${id}/${id}/outs/assembly \
--outprefix=${id}_pseudohap \
--minsize=10000 \
--headers=short


" >> $file_name

qsub $file_name
