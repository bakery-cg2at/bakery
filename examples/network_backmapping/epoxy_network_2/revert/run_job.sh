#!/bin/bash -l
#PBS -l mem=32gb
#PBS -l walltime=72:00:00
#PBS -o Output.job
#PBS -j oe
#PBS -l nodes=1:ppn=20
#PBS -M jakub.krajniak@cs.kuleuven.be
#PBS -A lp_polymer_goa_project

module purge

cd $PBS_O_WORKDIR

module load espressopp/adress
module load pyh5md
module load h5py/intel-gpfs
# Set up OpenMPI environment
n_proc=$(cat $PBS_NODEFILE | wc -l)
n_node=$(cat $PBS_NODEFILE | uniq | wc -l)

RES_STEP="`cat RES`"

EQ=5000
LONG=50000
TEMPERATURE=800.0
GAMMA=0.5
SKIN=0.1

date > ${OUTPUT_PREFIX}.log
echo "python -u start_sim.py --conf conf.gro --res_rate $RES_STEP --int_step 1000 --thermostat_gamma $GAMMA --long $LONG --eq $EQ --rng_seed ${SEED} --temperature $TEMPERATURE --output_prefix ${OUTPUT_PREFIX} --skin $SKIN --top topol.top  &>> ${OUTPUT_PREFIX}.log" >> ${OUTPUT_PREFIX}.log
python -u start_sim.py --conf conf.gro --res_rate $RES_STEP --temperature $TEMPERATURE --int_step 1000 --thermostat_gamma $GAMMA --long $LONG --eq $EQ --rng_seed ${SEED} --output_prefix ${OUTPUT_PREFIX} --skin $SKIN --top topol.top &>> ${OUTPUT_PREFIX}.log
date >> ${OUTPUT_PREFIX}.log
