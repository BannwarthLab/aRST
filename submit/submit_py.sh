#!/bin/bash

# ask for eight cores

#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=4
#SBATCH --time=48:00:00
#SBATCH --signal=B:USR1@120

#SBATCH --mem-per-cpu=1G   #M is the default and can therefore be omitted, but could also be K(ilo)|G(iga)|T(era)

# name the job
#SBATCH --job-name=react_search
# declare the merged STDOUT/STDERR file
#SBATCH --output=output.%J.txt

### beginning of executable commands
# Note: the OMP_NUM_THREADS envvar is set automatically - do not overwrite!
export MKL_NUM_THREADS=10
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

script_path="/home/guests/ychen/meow/autoreact"

# clean up old results
for file in $(ls $SLURM_SUBMIT_DIR);do
  if [[ -d $file ]];then
    rm -r $file
  elif [[ "$file" == *.json ]]; then
    rm $file
  fi
done

# upload script and input
for file in $(ls $script_path/*.py);do
  cp $file $tdir
done
for file in $(ls $SLURM_SUBMIT_DIR/*);do
  cp $file $tdir
done

cd $tdir

module load openmpi/4.1.1
module load orca/5.0.2
module load xtb/6.5.1
module load excstda/1.2
module load libraries/libgfortran/5


source "/home/guests/ychen/mb/env/bin/activate"

python3 main.py > py.out 2>&1 &
wait

shopt -s extglob
rm output*

for file in $(ls $tdir/*) ; do
  if [[ -d $file ]];then
    cp -r $file $SLURM_SUBMIT_DIR
  elif [[ "$file" == *.json ]]; then
    cp $file $SLURM_SUBMIT_DIR
  fi
done

rm -rf $tdir
cd  $SLURM_SUBMIT_DIR



      



