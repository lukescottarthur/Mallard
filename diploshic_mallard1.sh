#!/bin/bash
#SBATCH --job-name=diploshic_mallard1                      # Job name 
#SBATCH --partition=batch                           # Partition name 
#SBATCH --ntasks=1                                  # 
#SBATCH --cpus-per-task=16                           # CPU core count per task
#SBATCH --mem=32G                                    # Memory per node
#SBATCH --time=24:00:00                              # Time limit hrs:mins:secs
#SBATCH --output=/home/las80898/Mallard/%x_%j.out  ####change to your path
#SBATCH --error=/home/las80898/Mallard/%x_%j.error ####change to your path
#SBATCH --mail-user=las80898@uga.edu                #####change to your path
#SBATCH --mail-type=END,FAIL                        # Mail events 

# Load conda environment 
CONDA_BASE=$(conda info --base)
source ${CONDA_BASE}/etc/profile.d/conda.sh
conda activate diplo_env

# Navigate to diploSHIC directory (change to your directory)
cd /home/las80898/diploSHIC

# chr17 name
# chr17 length 215745

# make directory if necessary in diploshic (mine is 'mallardtest1')
# add vcf.gz file to your directory
# move discoal simulations into your directory

#calculate feature vectors from simulations
for f in mallardtest1/*_mallard1.gz; do diploSHIC fvecSim diploid $f $f.fvec --totalPhysLen 55000
done

# make balanced training set
mkdir rawFVFiles
mv mallardtest1/*.fvec rawFVFiles/
mkdir mallardTrainingSets
diploSHIC makeTrainingSets rawFVFiles/neutral_mallard1.gz.fvec rawFVFiles/soft rawFVFiles/hard 5 0,1,2,3,4,6,7,8,9,10 mallardTrainingSets/

# train!
diploSHIC train mallardTrainingSets/ mallardTrainingSets/ bfsModel

# feature vectors for real data
diploSHIC fvecVcf diploid mallardtest1/chr17merge.vcf.gz chr17 215745 mallardtest1/chr17merge.fvec --winSize 55000

# move output
mv mallardtest1/chr17merge.fvec rawFVFiles/chr17merge.fvec

# prediction on empirical feature vectors
diploSHIC predict bfsModel.json bfsModel.weights.h5 rawFVFiles/chr17merge.fvec mallardtest1_results