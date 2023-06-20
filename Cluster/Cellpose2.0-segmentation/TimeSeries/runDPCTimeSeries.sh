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
DIR=$(sed -n "${SLURM_ARRAY_TASK_ID}p" caseListTimeSeries)


# 1) Compute Node: Create virtual environment (node-local virtual environment) inside the job
module load gcc opencv python/3.10
virtualenv --no-download --clear $SLURM_TMPDIR/env && source $SLURM_TMPDIR/env/bin/activate
pip install cellpose

# 2) Compute Node: Run Cellpose2.0 on all images of the particular unique ID or UID (e.g. SCD-001, SCT-009, etc.) - based on DIR and the current array task ID
# The directory for UID is saved in variable $DIR, which is from the file caseListTimeSeries

# Loop through directories “CSID”  or coverslip IDs that start with Ai, Aii, Bi, Bii, and so on
# CSIDdir contains the full path of the CSID directory
for CSIDdir in /home/shrestha/scratch/ImagesCanada/DPC/TimeSeries/$DIR/[ABCDEFGNSTRabcdefgnstr]*/
do
        # Save the directory structure for UID/CSID in variable UtCSID
        UtCSID=${CSIDdir#*DPC/}
        #echo $UtCSID
        # Create directory (including intermediate directories - of UID, CSID - this is using ‘-p’) in OutlinesCanada folder if they do not exist
        # Directory OutlinesCanada/TimeSeries/UID/CSID is where the output files are saved, so any subdirectories not present are created
        mkdir -p /home/shrestha/projects/def-stoeber/shrestha/Cellpose/OutlinesCanada/$UtCSID
        # Run Cellpose on CSIDdir and save in savedir
        srun python -m cellpose --dir $CSIDdir --savedir /home/shrestha/projects/def-stoeber/shrestha/Cellpose/OutlinesCanada/$UtCSID --pretrained_model /home/shrestha/projects/def-stoeber/shrestha/Cellpose/models/sickle20230404  --diameter 0  --no_npy --save_txt --verbose --use_gpu
done
