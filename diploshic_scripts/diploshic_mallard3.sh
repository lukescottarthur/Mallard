#!/bin/bash
#SBATCH --job-name=diploshic_mallard3                      # Job name 
#SBATCH --partition=batch                           # Partition name 
#SBATCH --ntasks=1                                  # 
#SBATCH --cpus-per-task=16                           # CPU core count per task
#SBATCH --mem=32G                                    # Memory per node
#SBATCH --time=12:00:00                              # Time limit hrs:mins:secs
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

###### this is to try the pipeline with unmerged vcf files

# make new directory
# mkdir mallardtest3

# copy vcf files into directory
# cp /home/las80898/mallard_data/chr17.WMAIDreheader.vcf.gz mallardtest3/
# cp /home/las80898/mallard_data/chr17.GFMIDreheader.vcf.gz mallardtest3/

# feature vectors for real data WMA 
diploSHIC fvecVcf diploid mallardtest3/chr17.WMAIDreheader.vcf.gz chr17 215745 mallardtest3/chr17.WMA.fvec --winSize 55000

# feature vectors for real data GFM
diploSHIC fvecVcf diploid mallardtest3/chr17.GFMIDreheader.vcf.gz chr17 215745 mallardtest3/chr17.GFM.fvec --winSize 55000

# prediction on empirical feature vectors WMA
diploSHIC predict bfsModel.json bfsModel.weights.h5 mallardtest3/chr17.WMA.fvec mallardtest3WMA_results

# prediction on empirical feature vectors GFM
diploSHIC predict bfsModel.json bfsModel.weights.h5 mallardtest3/chr17.GFM.fvec mallardtest3GFM_results