#!/bin/bash

id=$1
time=$2
ppn=$3
kind=$4

echo "#!/bin/bash

#PBS -S /bin/bash
#PBS -A ACF-UTHSC0013
#PBS -l nodes=1:ppn=$ppn,walltime=$time
#PBS -m abe
#PBS -M hanyou@umich.edu
#PBS -j oe
#PBS -l $kind
#PBS -N $id

"
