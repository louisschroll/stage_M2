#!/bin/sh
#SBATCH --job-name=step1_model_selection
#SBATCH -N 1
#SBATCH -n 3
#SBATCH --mem-per-cpu=5G
#SBATCH --partition=cefe

### Email
#SBATCH --mail-user=louis.schroll@cefe.cnrs.fr         
#SBATCH --mail-type=ALL

### Output & Error
#SBATCH --output=/lustre/schrolll/3.results/1.test_quadratic_effect.out
#SBATCH --error=/lustre/schrolll/3.results/1.test_quadratic_effect.err              


echo « Running on: $SLURM_NODELIST »

srun /lustre/schrolll/2.code/pt1_spOccupancy/1.test_quadratic_effect.R