#!/bin/bash
#SBATCH --time=5:00:00
#SBATCH --mem=16GB
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --output=slurm/slurm_%j.out
#SBATCH --mail-type=END
#SBATCH --mail-user=azc211@nyu.edu

# Load modules
module purge
# Make sure bashrc is executed
source ~/.bashrc

# cd to repo directory
REPO=/home/azc211/github/ihdp
cd $REPO

# Command line arguments
export SCRIPT=$1
shift
export ARGS=$@

# Run program
echo "R $SCRIPT $ARGS"
singularity exec /beegfs/work/public/singularity/r-3.6.1.sif Rscript $SCRIPT $ARGS

