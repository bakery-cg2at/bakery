#!/bin/bash -le
#PBS -N em_annling
#PBS -l walltime=24:00:00
#PBS -o Output.job
#PBS -j oe
#PBS -l nodes=1:ppn=20
#PBS -M jakub.krajniak@cs.kuleuven.be
#PBS -A  lp_sim_interpoco

module purge
module load GROMACS

cd $PBS_O_WORKDIR

# Set up OpenMPI environment
n_proc=$(cat $PBS_NODEFILE | wc -l)
n_node=$(cat $PBS_NODEFILE | uniq | wc -l)

LOG="${PBS_O_WORKDIR}/${PBS_JOBID}.log"

MD="em.mdp nvt.mdp"
FIRST=1
LASTDIR=""

for mdp in $MD; do
    workdir=`basename $mdp .mdp`
    mkdir $workdir
    if [ "$FIRST" == "1" ]; then
        cp -v conf.gro $workdir
        FIRST=0
    else
        cp -v $LASTDIR/confout.gro $workdir/conf.gro
    fi
    cp $mdp $workdir
    cp topol.top $workdir
    cd $workdir
    gmx_mpi grompp -f $mdp -maxwarn 1 -v && gmx_mpi mdrun -ntomp $n_proc -v
    LASTDIR=$workdir
    cd ..
done
