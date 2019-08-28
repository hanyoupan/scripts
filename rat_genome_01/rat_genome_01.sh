#!/bin/bash

source /home/hanyou/script/envrionment/hanyous_modules_02.source

process=$1
runID=$2

home="/home/hanyou"
scratch="/scratch/junzli_flux/hanyou"


if [ "$process" == "picard_sortsam" ]; then

    time="100:00:00"
    memory="10gb"
    heapSize="9g"

    cd "$scratch/data2/bwa_sampe/$runID"
    for ff in `ls *.sam`; do

	bash "$home/script/pbs_header.sh" "`date +"%Y.%m.%d"`.$process.$runID.$ff.pbs" $time $memory > "$home/script/$runID/`date +"%Y.%m.%d"`.$process.$runID.$ff.pbs"
	echo "mkdir -p $scratch/data2/picard_sortsam/$runID
cd $scratch/data2/picard_sortsam/$runID
java -Xmx$heapSize -jar \$PICARD_JARS/SortSam.jar INPUT=$scratch/data2/bwa_sampe/$runID/$ff OUTPUT=$scratch/data2/picard_sortsam/$runID/$ff.bam SORT_ORDER=coordinate TMP_DIR=$scratch/tmp/ VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true
" >> "$home/script/$runID/`date +"%Y.%m.%d"`.$process"".$runID.$ff.pbs"
	cd "$home/script/$runID"
	qsub "$home/script/$runID/`date +"%Y.%m.%d"`.$process.$runID.$ff.pbs"

    done

fi


if [ "$process" == "picard_markduplicates" ]; then

    time="200:00:00"
    memory="20gb"
    heapSize="18g"

    cd "$scratch/data2/picard_sortsam/$runID"
    for ff in `ls *.bam`; do

        bash "$home/script/pbs_header.sh" "`date +"%Y.%m.%d"`.$process.$runID.$ff.pbs" $time $memory > "$home/script/$runID/`date +"%Y.%m.%d"`.$process.$runID.$ff.pbs"
        echo "mkdir -p $scratch/data2/picard_markduplicates/$runID
cd $scratch/data2/picard_markduplicates/$runID
java -Xmx$heapSize -jar \$PICARD_JARS/MarkDuplicates.jar \
INPUT=$scratch/data2/picard_sortsam/$runID/$ff \
OUTPUT=$scratch/data2/picard_markduplicates/$runID/$ff.bam \
METRICS_FILE=$scratch/data2/picard_markduplicates/$runID/$ff.txt \
ASSUME_SORTED=true \
OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 \
TMP_DIR=$scratch/tmp \
REMOVE_DUPLICATES=false \
VALIDATION_STRINGENCY=SILENT \
CREATE_INDEX=true
" >> "$home/script/$runID/`date +"%Y.%m.%d"`.$process"".$runID.$ff.pbs"
        cd "$home/script/$runID"
        qsub "$home/script/$runID/`date +"%Y.%m.%d"`.$process.$runID.$ff.pbs"

    done

fi


if [ "$process" == "picard_markduplicates-remove" ]; then

    time="200:00:00"
    memory="30gb"
    heapSize="27g"

    cd "$scratch/data2/picard_sortsam/$runID"
    for ff in `ls *.bam`; do

        bash "$home/script/pbs_header.sh" "`date +"%Y.%m.%d"`.$process.$runID.$ff.pbs" $time $memory > "$home/script/$runID/`date +"%Y.%m.%d"`.$process.$runID.$ff.pbs"
        echo "mkdir -p $scratch/data2/$process/$runID
cd $scratch/data2/$process/$runID
java -Xmx$heapSize -jar \$PICARD_JARS/MarkDuplicates.jar \\
INPUT=$scratch/data2/picard_sortsam/$runID/$ff \\
OUTPUT=$scratch/data2/$process/$runID/$ff.bam \\
METRICS_FILE=$scratch/data2/$process/$runID/$ff.txt \\
ASSUME_SORTED=true \\
OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 \\
TMP_DIR=$scratch/tmp \\
REMOVE_DUPLICATES=true \\
VALIDATION_STRINGENCY=SILENT \\
CREATE_INDEX=true
" >> "$home/script/$runID/`date +"%Y.%m.%d"`.$process"".$runID.$ff.pbs"
        cd "$home/script/$runID"
#        qsub "$home/script/$runID/`date +"%Y.%m.%d"`.$process.$runID.$ff.pbs"

    done

fi
