#!/bin/bash
for i in {21..40}
do
mkdir /lustre/haven/proj/UTHSC0013/UM/projects/rat_genome_03/long_ranger/20_hrlr_rn6/vcfs/"5739-JL-00"$i/
cp /lustre/haven/proj/UTHSC0013/UM/projects/rat_genome_03/long_ranger/20_hrlr_rn6/"5739-JL-00"$i/"5739-JL-00"$i/outs/phased_variants.vcf.gz /lustre/haven/proj/UTHSC0013/UM/projects/rat_genome_03/long_ranger/20_hrlr_rn6/vcfs/"5739-JL-00"$i/
cd /lustre/haven/proj/UTHSC0013/UM/projects/rat_genome_03/long_ranger/20_hrlr_rn6/vcfs/"5739-JL-00"$i/
gzip -d phased_variants.vcf.gz
mv phased_variants.vcf "5739-JL-00"$i".vcf"
done
