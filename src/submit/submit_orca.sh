#!/bin/bash

# ask for eight cores

#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=2
#SBATCH --time=12:00:00
#SBATCH --signal=B:USR1@120

#SBATCH --mem-per-cpu=1G   #M is the default and can therefore be omitted, but could also be K(ilo)|G(iga)|T(era)

# name the job
#SBATCH --job-name=orca
# declare the merged STDOUT/STDERR file
#SBATCH --output=output.%J.txt

### beginning of executable commands
# Note: the OMP_NUM_THREADS envvar is set automatically - do not overwrite!
#export MKL_NUM_THREADS=10
#export OMP_NUM_THREADS=10

your_cleanup_function() {
        echo "function your_cleanup_function called at $(date)"
        # do whatever cleanup you want here
        echo "Copying files to working directory, then cleaning up '$tdir'"
        rm $tdir/*tmp*
        mkdir -p $SLURM_SUBMIT_DIR/dump
        cp ./*  $SLURM_SUBMIT_DIR/dump && rm -rf $tdir
} # call your_cleanup_function once we receive USR1 signal

trap 'your_cleanup_function' USR1 SIGTERM SIGKILL

tdir=$(mktemp -d $TMPDIR/XXXXXXXXXXXXXXXXXX)

echo 'temporary directory were computations are run:' $tdir

echo $SLURM_SUBMIT_DIR
for f in $(ls $SLURM_SUBMIT_DIR/*) ; do
        cp $f $tdir
done

cd $tdir

module load openmpi/4.1.1
module load orca/5.0.2

### ORCA
$ORCA orca.inp  > orca.out && tac orca.out | grep -m 1 'FINAL SINGLE' |awk '{print $5}' > ENERGY &
wait

shopt -s extglob
rm output*
cp -r ./!(*tmp*)  $SLURM_SUBMIT_DIR/ && rm -rf $tdir

cd $SLURM_SUBMIT_DIR/
rm output*
