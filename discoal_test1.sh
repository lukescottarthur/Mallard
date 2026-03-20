#!/bin/bash
#SBATCH --job-name=discoal_test_1                      # Job name 
#SBATCH --partition=batch                           # Partition name 
#SBATCH --ntasks=1                                  # 
#SBATCH --cpus-per-task=16                           # CPU core count per task
#SBATCH --mem=32G                                    # Memory per node
#SBATCH --time=06:30:00                              # Time limit hrs:mins:secs
#SBATCH --output=/home/las80898/Mallard/%x_%j.out  ####change to your path
#SBATCH --error=/home/las80898/Mallard/%x_%j.error ####change to your path
#SBATCH --mail-user=las80898@uga.edu                #####change to your path
#SBATCH --mail-type=END,FAIL                        # Mail events 

cd /home/las80898/discoal

# number of chromosomes (sample size): 162
# number of replicates: 2000
# number of discrete sites in sequence: 55000
# -Pt (prior on theta) lower boundary: 20; upper boundary: 2000
# -ws (stochastic sweep): 0 generations ago
# -x (location of sweep on chromosome)
# -Pf (prior on standing variation) lower boundary: 0; upper boundary: 0.2


# 11x hard sweeps with different -x values
./discoal 162 2000 55000 -Pt 20 2000 -ws 0 -x 0.0454545 > mallardtest1/hard_0_mallard1.gz
./discoal 162 2000 55000 -Pt 20 2000 -ws 0 -x 0.1363636 > mallardtest1/hard_1_mallard1.gz
./discoal 162 2000 55000 -Pt 20 2000 -ws 0 -x 0.2272727 > mallardtest1/hard_2_mallard1.gz
./discoal 162 2000 55000 -Pt 20 2000 -ws 0 -x 0.3181818 > mallardtest1/hard_3_mallard1.gz
./discoal 162 2000 55000 -Pt 20 2000 -ws 0 -x 0.4090909 > mallardtest1/hard_4_mallard1.gz
./discoal 162 2000 55000 -Pt 20 2000 -ws 0 -x 0.5 > mallardtest1/hard_5_mallard1.gz
./discoal 162 2000 55000 -Pt 20 2000 -ws 0 -x 0.5909090 > mallardtest1/hard_6_mallard1.gz
./discoal 162 2000 55000 -Pt 20 2000 -ws 0 -x 0.6818181 > mallardtest1/hard_7_mallard1.gz
./discoal 162 2000 55000 -Pt 20 2000 -ws 0 -x 0.7727272 > mallardtest1/hard_8_mallard1.gz
./discoal 162 2000 55000 -Pt 20 2000 -ws 0 -x 0.8636363 > mallardtest1/hard_9_mallard1.gz
./discoal 162 2000 55000 -Pt 20 2000 -ws 0 -x 0.9545454 > mallardtest1/hard_10_mallard1.gz

# 11x soft sweeps with different -x values
./discoal 162 2000 55000 -Pt 20 2000 -ws 0 -Pf 0 0.2 -x 0.0454545 > mallardtest1/soft_0_mallard1.gz
./discoal 162 2000 55000 -Pt 20 2000 -ws 0 -Pf 0 0.2 -x 0.1363636 > mallardtest1/soft_1_mallard1.gz
./discoal 162 2000 55000 -Pt 20 2000 -ws 0 -Pf 0 0.2 -x 0.2272727 > mallardtest1/soft_2_mallard1.gz
./discoal 162 2000 55000 -Pt 20 2000 -ws 0 -Pf 0 0.2 -x 0.3181818 > mallardtest1/soft_3_mallard1.gz
./discoal 162 2000 55000 -Pt 20 2000 -ws 0 -Pf 0 0.2 -x 0.4090909 > mallardtest1/soft_4_mallard1.gz
./discoal 162 2000 55000 -Pt 20 2000 -ws 0 -Pf 0 0.2 -x 0.5 > mallardtest1/soft_5_mallard1.gz
./discoal 162 2000 55000 -Pt 20 2000 -ws 0 -Pf 0 0.2 -x 0.5909090 > mallardtest1/soft_6_mallard1.gz
./discoal 162 2000 55000 -Pt 20 2000 -ws 0 -Pf 0 0.2 -x 0.6818181 > mallardtest1/soft_7_mallard1.gz
./discoal 162 2000 55000 -Pt 20 2000 -ws 0 -Pf 0 0.2 -x 0.7727272 > mallardtest1/soft_8_mallard1.gz
./discoal 162 2000 55000 -Pt 20 2000 -ws 0 -Pf 0 0.2 -x 0.8636363 > mallardtest1/soft_9_mallard1.gz
./discoal 162 2000 55000 -Pt 20 2000 -ws 0 -Pf 0 0.2 -x 0.9545454 > mallardtest1/soft_10_mallard1.gz

# neutral sweep
./discoal 162 2000 55000 -Pt 20 2000 > mallardtest1/neutral_mallard1.gz