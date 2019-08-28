#!/bin/bash

# Environment: source hanyous_modules_02.sh

process=$1
runID=$2

home="/home/hanyou"
scratch="/scratch/junzli_flux/hanyou"

refDir="$scratch/resource/reference/rat_genome_02"
ref="$scratch/resource/reference/rat_genome_02/rn6.clean.fa"


if [ "$process" == "long_ranger_mkref" ]; then

    time="3:00:00"
    memory="11gb"

    bash "$home/script/pbs_header.sh" "`date +"%Y.%m.%d"`.$process.$runID.pbs" $time $memory > "$home/script/$runID/`date +"%Y.%m.%d"`.$process.$runID.pbs"
    echo "cd $refDir
longranger mkref $ref
" >> "$home/script/$runID/`date +"%Y.%m.%d"`.$process"".$runID.pbs"

    cd "$home/script/$runID"
    qsub "$home/script/$runID/`date +"%Y.%m.%d"`.$process.$runID.pbs"

fi


