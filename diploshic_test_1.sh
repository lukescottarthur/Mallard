#!/bin/bash
#SBATCH --job-name=diploshic_test_1                      # Job name 
#SBATCH --partition=batch                           # Partition name 
#SBATCH --ntasks=1                                  # 
#SBATCH --cpus-per-task=8                           # CPU core count per task
#SBATCH --mem=16G                                    # Memory per node
#SBATCH --time=03:30:00                              # Time limit hrs:mins:secs
#SBATCH --output=/home/las80898/diploSHIC/%x_%j.out  ####change to your path
#SBATCH --error=/home/las80898/diploSHIC/%x_%j.error ####change to your path
#SBATCH --mail-user=las80898@uga.edu                # Where to send mail
#SBATCH --mail-type=BEGIN,END,FAIL                        # Mail events 

# Load conda environment 
CONDA_BASE=$(conda info --base)
source ${CONDA_BASE}/etc/profile.d/conda.sh
conda activate diplo_env

# Navigate to diploSHIC directory (change to your directory)
cd /home/las80898/diploSHIC

for f in exampleApplication/*.msOut.gz; do diploSHIC fvecSim diploid $f $f.diploid.fvec --totalPhysLen 55000 --maskFileName exampleApplication/Anopheles-gambiae-PEST_CHROMOSOMES_AgamP3.accessible.fa.gz --chrArmsForMasking 3R 
done

# make balanced training set
mkdir rawFVFiles
mv exampleApplication/*.fvec rawFVFiles/
mkdir trainingSets
diploSHIC makeTrainingSets rawFVFiles/neut.msOut.gz.diploid.fvec rawFVFiles/soft rawFVFiles/hard 5 0,1,2,3,4,6,7,8,9,10 trainingSets/

# train! (kind of long)
diploSHIC train trainingSets/ trainingSets/ bfsModel

# feature vectors for real data
diploSHIC fvecVcf diploid exampleApplication/ag1000g.phase1.ar3.pass.biallelic.3R.vcf.28000000-29000000.gz 3R 53200684 exampleApplication/ag1000g.phase1.ar3.pass.biallelic.3R.vcf.28000000-29000000.gz.diploid.fvec --targetPop BFS --sampleToPopFileName exampleApplication/samples_pops.txt --winSize 55000 --maskFileName exampleApplication/Anopheles-gambiae-PEST_CHROMOSOMES_AgamP3.accessible.fa.gz

# move output
mv exampleApplication/ag1000g.phase1.ar3.pass.biallelic.3R.vcf.28000000-29000000.gz.diploid.fvec rawFVFiles/ag1000g.phase1.ar3.pass.biallelic.3R.vcf.28000000-29000000.gz.diploid.fvec

# prediction on empirical feature vectors
diploSHIC predict bfsModel.json bfsModel.weights.h5 rawFVFiles/ag1000g.phase1.ar3.pass.biallelic.3R.vcf.28000000-29000000.gz.diploid.fvec mossie.preds