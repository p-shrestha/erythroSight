#!/bin/bash
#SBATCH --account=def-stoeber
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --gpus-per-node=v100l:1
#SBATCH --mem=8G
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
"ABeta-002"
"SS-011"
"AA-021"
"AA-022"
"AA-024"
"SS-012"
"AS-016"
"AA-012"
"AA-011"
"AS-034"
"AS-027"
"AS-031"
"AS-011"
"AA-004"
"SS-017"
"AS-015"
"AS-017"
"AA-015"
"ABeta-010"
"AA-010"
"SBeta-007"
"ABeta-008"
"AS-009"
"AA-014"
"ABeta-020"
"ABeta-014"
"AA-016"
"SBeta-009"
"AA-018"
"SS-009"
"SS-010"
"AA-017"
"AS-020"
"SS-013"
"AS-030"
"AS-012"
"AA-013"
"SS-018"
"AS-036"
"AS-028"
"AS-021"
"AS-022"
"AS-013"
"AS-025"
"ABeta-015"
"ABeta-022"
"ABeta-011"
"ABeta-009"
"AA-020"
"ABeta-021"
"AA-019"
"ABeta-016"
"AS-026"
"ABeta-024"
"ABeta-023"
"AS-014"
"AA-023"
"AS-029"
"ABeta-017"
"AS-018"
"SS-015"
"SBeta-002"
"AA-003"
"ABeta-007"
"AS-024"
"AS-033"
"AA-001"
"AA-002"
"AS-032"
"SS-008"
"AA-009"
"AS-019"
"AS-010"
"SS-002"
"SS-006"
"AA-008"
"AS-004"
"AS-003"
"SS-016"
"ABeta-001"
"AS-002"
"SBeta-003"
"ABeta-006"
"SS-005"
"AS-005"
"AA-006"
"SBeta-001"
"ABeta-003"
"SBeta-004"
"AA-005"
"ABeta-004"
"ABeta-005"
"SS-001"
"SBeta-005"
"AA-007"
"SS-007"
"ABeta-013"
"SBeta-006"
"AS-035"
"SS-014"
"SBeta-010"
"ABeta-019"
"SS-003"
"AS-007"
"ABeta-012"
"AS-006"
"AS-008"
"SS-004"
"SBeta-008"
"ABeta-018"
"AS-023"
"AS-001"
"Test"
"Test1"
)

# 4) Compute Node (virtual): Run ImageJ (headless mode)

# Check if SLURM_ARRAY_TASK_ID is within the valid range
if (( $SLURM_ARRAY_TASK_ID >= 1 && $SLURM_ARRAY_TASK_ID <= 114 )); then
  # Execute the command corresponding to the SLURM_ARRAY_TASK_ID
  command="ImageJ-linux64 --headless --console --run /home/shrestha/projects/def-stoeber/shrestha/MorphologicalParameters/Scripts/imagej_roi_converter_All_CC.py 'UID=\"${uids[$SLURM_ARRAY_TASK_ID - 1]}\" ,dataPath=\"/home/shrestha/projects/def-stoeber/shrestha/Cellpose/Outlines/\" , resultsPath=\"/home/shrestha/projects/def-stoeber/shrestha/MorphologicalParameters/Nepal/\"'"
  eval "$command"
else
  echo "Invalid case: $SLURM_ARRAY_TASK_ID"
fi
