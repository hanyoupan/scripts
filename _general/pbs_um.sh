$#!/bin/bash

jobID=$1
time=$2
memory=$3
ppn=$4

if [ "$ppn" == "" ]; then
    ppn=1
fi

echo "#!/bin/bash

#This tells the scheduler what account to use.  Do not change.
#PBS -A junzli_flux
#This is the job name.  Feel free to change this at will to whatever you need.
#PBS -N $jobID
#This denotes the queue that the job should be run in.  Do not change if you want to use the flux dedicated nodes
#PBS -q flux
#The next two denotes the email address to send jobs to, and under what conditions to send that email.
#PBS -M hanyou@umich.edu
#This line says (a) send email if the job fails, (b) when the job starts, and (e) when the job ends.
#PBS -m abe
#This line combine e file into o file.
#PBS -j oe
#This line sends all environment variables on the login node.
#PBS -V
#This denotes the number of nodes and processors that the job should be run on. The max for ppn is currently 12.
#Walltime is denoted by hh:mm:ss, and hours can be no more than 376 (or 28 days).
#PBS -l nodes=1:ppn=$ppn,pmem=$memory,walltime=$time
#PBS -l qos=flux
"

