#!/bin/bash
#SBATCH --job-name=diploshic_test_1                      # Job name 
#SBATCH --partition=batch                           # Partition name 
#SBATCH --ntasks=1                                  # 1 task (process)
#SBATCH --cpus-per-task=8                           # CPU core count per task
#SBATCH --mem=16G                                    # Memory per node
#SBATCH --time=08:00:00                              # Time limit hrs:mins:secs
#SBATCH --output=/home/las80898/diploSHIC/%x_%j.out  
#SBATCH --error=/home/las80898/diploSHIC/%x_%j.error 
#SBATCH --mail-user=las80898@uga.edu                # Where to send mail
#SBATCH --mail-type=BEGIN,END,FAIL                        # Mail events 

# Load conda environment 
CONDA_BASE=$(conda info --base)
source ${CONDA_BASE}/etc/profile.d/conda.sh
conda activate diplo_env

# Change directory 
cd /home/las80898/diploSHIC

for f in exampleApplication/*.msOut.gz 
do 
diploSHIC fvecSim diploid $f $f.diploid.fvec --totalPhysLen 55000 --maskFileName exampleApplication/Anopheles-gambiae-PEST_CHROMOSOMES_AgamP3.accessible.fa.gz --chrArmsForMasking 3R 
& 
done
