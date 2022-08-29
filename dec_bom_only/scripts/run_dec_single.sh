#!/bin/bash
#SBATCH --partition=shared
#SBATCH --threads-per-core=1
#SBATCH --mem=1024M
#SBATCH --mail-user=ctribble09@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --time=72:00:00
#SBATCH --ntasks=1

module load devel/Boost
cd dec_bom_only

../rb scripts/run_dec_single.Rev