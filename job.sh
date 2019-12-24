#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=5
#SBATCH --job-name=test
#SBATCH --output=test.out
#SBATCH --mem=2000000
#SBATCH --partition=himem

module load R/3.4.2-foss-2016a-X11-20160819
Rscript test_code.R