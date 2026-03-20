#!/bin/bash
#SBATCH --job-name=diploshic_mallard2                      # Job name 
#SBATCH --partition=batch                           # Partition name 
#SBATCH --ntasks=1                                  # 
#SBATCH --cpus-per-task=16                           # CPU core count per task
#SBATCH --mem=32G                                    # Memory per node
#SBATCH --time=4:00:00                              # Time limit hrs:mins:secs
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


# feature vectors for real data WMA 
diploSHIC fvecVcf diploid mallardtest1/chr17merge.vcf.gz chr17 215745 mallardtest1/chr17mergeWMA.fvec --targetPop WMA --sampleToPopFileName mallardtest1/mallard_pop_id.txt --winSize 55000

# move output
mv mallardtest1/chr17mergeWMA.fvec rawFVFiles/chr17mergeWMA.fvec

# prediction on empirical feature vectors
diploSHIC predict bfsModel.json bfsModel.weights.h5 rawFVFiles/chr17mergeWMA.fvec mallardtest1WMA_results

# feature vectors for real data GFM
diploSHIC fvecVcf diploid mallardtest1/chr17merge.vcf.gz chr17 215745 mallardtest1/chr17mergeGFM.fvec --targetPop GFM --sampleToPopFileName mallardtest1/mallard_pop_id.txt --winSize 55000

# move output
mv mallardtest1/chr17mergeGFM.fvec rawFVFiles/chr17mergeGFM.fvec

# prediction on empirical feature vectors
diploSHIC predict bfsModel.json bfsModel.weights.h5 rawFVFiles/chr17mergeGFM.fvec mallardtest1GFM_results