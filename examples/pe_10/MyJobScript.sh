#!/bin/bash -l
#PBS -N R_0_100  
#PBS -l mem=8gb
#PBS -l walltime=20:00:00
#PBS -o Output.job
#PBS -j oe
#PBS -l nodes=1:ppn=20
module load espressopp/adress_new 
module load bakery-github-dev

cd $PBS_O_WORKDIR

mpirun -np 18 start_backmapping.py @params 
