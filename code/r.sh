#!/bin/bash
#SBATCH --time=4:00:00
#SBATCH --mem=4GB
#SBATCH --nodes=10
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1

# Load modules
module purge
module load r/intel/3.6.0

# cd to repo directory
REPO=/home/azc211/github/ihdp
cd $REPO

# Command line arguments
export R_SCRIPT=$1
shift
export ARGS=$@

# Run program
echo "Rscript $R_SCRIPT $ARGS"
Rscript $R_SCRIPT $ARGS
