#!/bin/bash
#SBATCH --account=def-stoeber
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=0
#SBATCH --mail-type=begin         # send mail when job begins
#SBATCH --mail-type=end           # send mail when job ends
#SBATCH --mail-type=fail          # send mail if job fails
#SBATCH --mail-type=array_tasks   # send mail for all job array tasks
#SBATCH --mail-user=pranavsh@mail.ubc.ca

echo "Starting task $SLURM_ARRAY_TASK_ID"

# 1) Compute Node: Create virtual environment (node-local virtual environment) inside the job
module load gcc opencv python/3.10
virtualenv --no-download --clear $SLURM_TMPDIR/env && source $SLURM_TMPDIR/env/bin/activate

# 2) Compute Node (virtual): Activate matlab
module load matlab/2023a

# 3) Define array of group names
gNames=(
"AA"
"ABeta"
"AS"
"SCD"
)

# 4) Compute Node (virtual): Run MATLAB program
matlab -nodisplay -r "groupName = \"${gNames[$SLURM_ARRAY_TASK_ID - 1]}\";FrequencyDistribution1D_100b40mp"
