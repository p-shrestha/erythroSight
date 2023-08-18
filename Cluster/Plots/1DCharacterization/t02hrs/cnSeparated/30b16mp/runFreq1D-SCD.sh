#!/bin/bash
#SBATCH --account=def-stoeber
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=0
#SBATCH --time=02:00:00
#SBATCH --mail-type=begin         # send mail when job begins
#SBATCH --mail-type=end           # send mail when job ends
#SBATCH --mail-type=fail          # send mail if job fails
#SBATCH --mail-user=pranavsh@mail.ubc.ca

# 1) Compute Node: Create virtual environment (node-local virtual environment) inside the job
module load gcc opencv python/3.10
virtualenv --no-download --clear $SLURM_TMPDIR/env && source $SLURM_TMPDIR/env/bin/activate

# 2) Compute Node (virtual): Activate matlab
module load matlab

# 3) Compute Node (virtual): Run MATLAB program
#matlab -nodisplay -r "groupName = \"AA\";tName = \"t02hrs\";FrequencyDistribution1D_30b16mp"
#matlab -nodisplay -r "groupName = \"ABeta\";tName = \"t02hrs\";FrequencyDistribution1D_30b16mp"
#matlab -nodisplay -r "groupName = \"AS\";tName = \"t02hrs\";FrequencyDistribution1D_30b16mp"
#matlab -nodisplay -r "groupName = \"SBeta\";tName = \"t02hrs\";FrequencyDistribution1D_30b16mp"
#matlab -nodisplay -r "groupName = \"SS\";tName = \"t02hrs\";FrequencyDistribution1D_30b16mp"
#matlab -nodisplay -r "groupName = \"N\";tName = \"t02hrs\";FrequencyDistribution1D_30b16mp"
#matlab -nodisplay -r "groupName = \"SCT\";tName = \"t02hrs\";FrequencyDistribution1D_30b16mp"
matlab -nodisplay -r "groupName = \"SCD\";tName = \"t02hrs\";FrequencyDistribution1D_30b16mp"
