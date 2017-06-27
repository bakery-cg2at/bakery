#!/bin/bash -le
#PBS -N em_annling
#PBS -l walltime=24:00:00
#PBS -o Output.job
#PBS -j oe
#PBS -l nodes=1:ppn=20
#PBS -M jakub.krajniak@cs.kuleuven.be
#PBS -A  lp_sim_interpoco

source switch_to_2015a

module purge
module load 2015a/espressopp/adress
module load bakery-github-dev

cd $PBS_O_WORKDIR

# Set up OpenMPI environment
n_proc=$(cat $PBS_NODEFILE | wc -l)
n_node=$(cat $PBS_NODEFILE | uniq | wc -l)
#mpdboot -f $PBS_NODEFILE -n $n_node -r ssh -v

LOG="${PBS_O_WORKDIR}/${PBS_JOBID}.log"

for alpha in 0.001 0.0001 0.00001 0.000001; do
    for s in {1..3}; do
        mkdir -p data_s${s}
        #mkdir -p gmx_${s}/${alpha}
        rng_seed=`shuf -n1 -i1000-10000`
        mpirun -n 18 start_backmapping.py @params --alpha $alpha --rng_seed $rng_seed --output_prefix data_s${s}/sim0 --trj_collect 10000 &> ${LOG}_${alpha}_$s
    done
done