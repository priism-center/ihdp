#!/bin/bash
#SBATCH --time=5:00:00
#SBATCH --mem=16GB
#SBATCH --nodes=10
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --output=slurm/slurm_%j.out
#SBATCH --mail-type=END
#SBATCH --mail-user=azc211@nyu.edu

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