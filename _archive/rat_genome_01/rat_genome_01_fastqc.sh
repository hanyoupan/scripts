#!/bin/bash

# Environment: source Bilges_modules.sh

process=$1
runID=$2

home="/home/hanyou"
scratch="/scratch/junzli_flux/hanyou"

if [ "$process" == "fastqc" ]; then

    time="5:00:00"
    memory="8gb"

    sequenceDir=$scratch"/data2/sequence/"$runID
    cd $sequenceDir

    for ff in `ls`; do

	bash $home"/script/PBSHeader.sh" "`date +"%Y.%m.%d"`"".$process""."$ff".pbs" $time $memory > "$home/script/$runID/""`date +"%Y.%m.%d"`"".$process""."$ff".pbs"
	mkdir -p "$scratch/data2/fastqc/$runID"
	echo "fastqc $sequenceDir/$ff --outdir $scratch/data2/fastqc/$runID/" >> "$home/script/$runID/""`date +"%Y.%m.%d"`"".$process""."$ff".pbs"
	
	cd "$home/script/$runID"
	qsub "$home/script/$runID/""`date +"%Y.%m.%d"`"".$process""."$ff".pbs"

	done

fi
