#!/bin/bash
#SBATCH --account=def-stoeber
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --gpus-per-node=v100l:1
#SBATCH --mem=100G
#SBATCH --mail-type=begin         # send mail when job begins
#SBATCH --mail-type=end           # send mail when job ends
#SBATCH --mail-type=fail          # send mail if job fails
#SBATCH --mail-type=array_tasks   # send mail for all job array tasks
#SBATCH --mail-user=pranavsh@mail.ubc.ca

echo "Starting task $SLURM_ARRAY_TASK_ID"

# 1) Compute Node: Create virtual environment (node-local virtual environment) inside the job
module load gcc opencv python/3.10
virtualenv --no-download --clear $SLURM_TMPDIR/env && source $SLURM_TMPDIR/env/bin/activate

# 2) Compute Node (virtual): Activate ImageJ
module load StdEnv/2020 fiji/20201104-1356

# 3) Define array of UIDs; Change this for different cases (e.g. for Nepal vs Canada data)
uids=(
"SCT-007"
"SCD-009"
"SCD-007"
"SCD-010"
"SCD-005"
"SCT-006"
"SCD-012"
"SCT-003"
"SCD-006"
"SCD-011"
"SCT-009"
"SCD-004"
"SCD-008"
"SCD-001"
"SCT-001"
)

# 4) Compute node: Copy python program in $SLURM_TMPDIR directory
cd $SLURM_TMPDIR
cp /home/shrestha/scratch/ScriptsCCScratch/TimeSeriesTIFF/createTIFFstacks.py .


# 5) Compute Node (virtual): Run ImageJ (headless mode)

# Check if SLURM_ARRAY_TASK_ID is within the valid range
if (( $SLURM_ARRAY_TASK_ID >= 1 && $SLURM_ARRAY_TASK_ID <= 30 )); then
  # Execute the command corresponding to the SLURM_ARRAY_TASK_ID
  command="ImageJ-linux64 --headless --console --run $SLURM_TMPDIR/createTIFFstacks.py 'UID=\"${uids[$SLURM_ARRAY_TASK_ID - 1]}\" ,dataPath=\"/home/shrestha/scratch/ImagesCanada/DPC/TimeSeries/\" , resultsPath=\"/home/shrestha/scratch/ImagesCanada/DPC/TimeSeries/TIFFstacks/\"'"
  eval "$command"
else
  echo "Invalid case: $SLURM_ARRAY_TASK_ID"
fi
