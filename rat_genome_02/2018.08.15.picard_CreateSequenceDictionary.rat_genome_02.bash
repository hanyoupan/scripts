#!/bin/bash

# Environment: source hanyous_modules_01.sh

cd /scratch/junzli_flux/hanyou/resource/reference/rat_genome_02
java -jar $PICARD_JARS/CreateSequenceDictionary.jar R=./refdata-rn6.clean/fasta/genome.fa O=./refdata-rn6.clean/fasta/genome.dict
