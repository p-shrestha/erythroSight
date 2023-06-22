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
"RS-N-005"
"SCD-012"
"TS-SCD-004"
"SCD-001"
"SCT-009"
"SCD-009"
"SCD-010"
"N-006"
"N-004"
"N-005"
"SCD-002"
"SCT-008"
"SCD-004"
"N-003"
"N-002"
"SCD-003"
"N-001"
"SCT-003"
"SCD-007"
"SCT-005"
"SCD-011"
"SCT-007"
"SCT-002"
"SCT-004"
"SCD-008"
"SCD-006"
"SCT-001"
"SCT-006"
"SCD-005"
)

# 4) Compute node: Make directories in the compute mode $SLURM_TMPDIR directory; and copy python program in current directory
cd $SLURM_TMPDIR
mkdir Outlines
mkdir MorphParam
cp /home/shrestha/projects/def-stoeber/shrestha/MorphologicalParameters/Scripts/imagej_roi_converter_CC_ArchiveCanadaR.py .

# 5) Login node and Compute Node (virtual): Create archive of Outlines and save it in compute node $SLURM_TMPDIR
cd /home/shrestha/projects/def-stoeber/shrestha/Cellpose/OutlinesCanada/
tar -cf $SLURM_TMPDIR/Outlines/${uids[$SLURM_ARRAY_TASK_ID - 1]}.tar ${uids[$SLURM_ARRAY_TASK_ID - 1]}

# Compute node: Extract .tar file in local storage
cd $SLURM_TMPDIR/Outlines
tar -xf $SLURM_TMPDIR/Outlines/${uids[$SLURM_ARRAY_TASK_ID - 1]}.tar

# 6) Compute Node (virtual): Run ImageJ (headless mode)

# Check if SLURM_ARRAY_TASK_ID is within the valid range
if (( $SLURM_ARRAY_TASK_ID >= 1 && $SLURM_ARRAY_TASK_ID <= 30 )); then
  # Execute the command corresponding to the SLURM_ARRAY_TASK_ID
  command="ImageJ-linux64 --headless --console --run $SLURM_TMPDIR/imagej_roi_converter_CC_ArchiveCanadaR.py 'UID=\"${uids[$SLURM_ARRAY_TASK_ID - 1]}\" ,dataPath=\"$SLURM_TMPDIR/Outlines/\" , resultsPath=\"$SLURM_TMPDIR/MorphParam/\"'"
  eval "$command"
else
  echo "Invalid case: $SLURM_ARRAY_TASK_ID"
fi

# The computations are done, so clean up the data set...
# Login node: Save the processed morphological parameters (archived) in the scratch space
cd $SLURM_TMPDIR/MorphParam/
tar -cf /home/shrestha/projects/def-stoeber/shrestha/MorphologicalParameters/Canada/${uids[$SLURM_ARRAY_TASK_ID - 1]}_MorphologicalParameters.tar ${uids[$SLURM_ARRAY_TASK_ID - 1]}
