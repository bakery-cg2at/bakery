#!/bin/bash -l
#PBS -N PE_423_NVT_1
#PBS -l walltime=24:00:00
#PBS -o Output.job
#PBS -j oe
#PBS -l nodes=1:ppn=24
#PBS -M jakub.krajniak@cs.kuleuven.be
#PBS -A lp_polymer_goa_project


module purge
module load GROMACS/5.0.4-intel-2014a-hybrid
module remove impi
module load impi
MDRUN="mdrun_mpi -v" 
GROMPP="grompp_mpi"

cd $PBS_O_WORKDIR

# Set up OpenMPI environment
n_proc=$(cat $PBS_NODEFILE | wc -l)
n_node=$(cat $PBS_NODEFILE | uniq | wc -l)
#mpdboot -f $PBS_NODEFILE -n $n_node -r ssh -v

mpirun -np $n_proc $MDRUN 

mpdallexit
