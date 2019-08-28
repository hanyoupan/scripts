#!/bin/bash

source /home/hanyou/script/envrionment/hanyous_modules_02.source

process=$1
runID=$2

home="/home/hanyou"
scratch="/scratch/junzli_flux/hanyou"

refDir="$scratch/resource/reference/rat_genome_01"
ref="$scratch/resource/reference/rat_genome_01/rn6.fa"


if [ "$process" == "bwa_index" ]; then

    time="10:00:00"
    memory="10gb"

    bash "$home/script/pbs_header.sh" "`date +"%Y.%m.%d"`.$process.$runID.pbs" $time $memory > "$home/script/$runID/`date +"%Y.%m.%d"`.$process.$runID.pbs"
    echo "cd $refDir
bwa index -a bwtsw $ref
" >> "$home/script/$runID/`date +"%Y.%m.%d"`.$process"".$runID.pbs"

    cd "$home/script/$runID"
    qsub "$home/script/$runID/`date +"%Y.%m.%d"`.$process.$runID.pbs"

fi


if [ "$process" == "bwa_aln" ]; then

    time="150:00:00"
    memory="4gb"

    cd "$scratch/data2/sequence/$runID"
    for ff in `ls *.gz`; do

	bash "$home/script/pbs_header.sh" "`date +"%Y.%m.%d"`.$process.$runID.$ff.pbs" $time $memory > "$home/script/$runID/`date +"%Y.%m.%d"`.$process.$runID.$ff.pbs"
	echo "cd $scratch/data2/sequence/$runID
mkdir -p $scratch/data2/bwa_aln/$runID
bwa aln -q 15 $ref $ff > $scratch/data2/bwa_aln/$runID/$ff.sai
" >> "$home/script/$runID/`date +"%Y.%m.%d"`.$process"".$runID.$ff.pbs"

	cd "$home/script/$runID"
	qsub "$home/script/$runID/`date +"%Y.%m.%d"`.$process.$runID.$ff.pbs"

    done
    
fi


if [ "$process" == "bwa_sampe" ]; then

    time="200:00:00"
    memory="10gb"

    sai1="/scratch/junzli_flux/hanyou/data2/bwa_aln/rat_genome_01/s-BN-N_HL53WCCXX_L2_1.clean.fq.gz.sai"
    sai2="/scratch/junzli_flux/hanyou/data2/bwa_aln/rat_genome_01/s-BN-N_HL53WCCXX_L2_2.clean.fq.gz.sai"
    fq1="/scratch/junzli_flux/hanyou/data2/sequence/rat_genome_01/s-BN-N_HL53WCCXX_L2_1.clean.fq.gz"
    fq2="/scratch/junzli_flux/hanyou/data2/sequence/rat_genome_01/s-BN-N_HL53WCCXX_L2_2.clean.fq.gz"
    samout="/scratch/junzli_flux/hanyou/data2/bwa_sampe/rat_genome_01/s-BN-N_HL53WCCXX_L2.clean.fq.gz.sam"

    bash "$home/script/pbs_header.sh" "`date +"%Y.%m.%d"`.$process.$runID.pbs" $time $memory > "$home/script/$runID/`date +"%Y.%m.%d"`.$process.$runID.pbs"
    echo "mkdir -p $scratch/data2/bwa_sampe/$runID
bwa sampe -r \"@RG\tID:Run_Novogene-Rats_BN-N\tSM:BN-N\tPL:ILLUMINA\tLB:BN-N\" $ref $sai1 $sai2 $fq1 $fq2 > $samout
" >> "$home/script/$runID/`date +"%Y.%m.%d"`.$process"".$runID.pbs"

    cd "$home/script/$runID"
    qsub "$home/script/$runID/`date +"%Y.%m.%d"`.$process.$runID.pbs"

fi
