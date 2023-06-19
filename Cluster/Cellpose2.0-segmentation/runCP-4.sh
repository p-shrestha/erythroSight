#!/bin/bash
#SBATCH --account=def-stoeber
#SBATCH --array=108-112%3         # Run array of 3 simultaneous jobs from 108 to 112 in case_list
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --gpus-per-node=v100l:1
#SBATCH --mem=8G
#SBATCH --time=5-05:00:00         # Time is set to 125 hours (5 days, 5 hours) for each job in the array
#SBATCH --mail-type=begin         # send mail when job begins
#SBATCH --mail-type=end           # send mail when job ends
#SBATCH --mail-type=fail          # send mail if job fails
#SBATCH --mail-type=array_tasks   # send mail for all job array tasks
#SBATCH --mail-user=pranavsh@mail.ubc.ca

echo "Starting task $SLURM_ARRAY_TASK_ID"
DIR=$(sed -n "${SLURM_ARRAY_TASK_ID}p" case_list)
#cd $DIR


# 1) Compute Node: Create virtual environment (node-local virtual environment) inside the job
module load gcc opencv python/3.10
virtualenv --no-download --clear $SLURM_TMPDIR/env && source $SLURM_TMPDIR/env/bin/activate
pip install cellpose

# 2) Compute Node: Run Cellpose2.0 on all images of the particular unique ID or UID (e.g. AA-001, AS-023, etc.) - based on DIR and the current array task ID
# The directory for UID is saved in variable $DIR, which is from the file case_list
# Loop through sub-directories “tname” (e.g. t0, t02hrs, etc.)
# tnamedir contains the full path of the tname directory
for tnamedir in ~/scratch/ImagesNepal/DPC/$DIR/t*/
do
        # Loop through sub-directories “CSID”  or coverslip IDs that start with Ai, Aii, Bi, Bii, and so on
        # CSIDdir contains the full path of the CSID directory
        for CSIDdir in $tnamedir[ABCDEFNSabcdefns]*/
        do
                # Save the directory structure for UID/tname/CSID in variable UtCSID
                UtCSID=${CSIDdir#*DPC/}
                #echo $UtCSID
                # Create directory (including intermediate directories - of UID, tname, CSID - this is using ‘-p’) in Outlines folder if they do not exist
                # Directory Outlines/UID/tname/CSID is where the output files are saved, so any subdirectories not present are created
                mkdir -p ~/projects/def-stoeber/shrestha/Cellpose/Outlines/$UtCSID
                # Run Cellpose on CSIDdir and save in savedir
                srun python -m cellpose --dir $CSIDdir --savedir ~/projects/def-stoeber/shrestha/Cellpose/Outlines/$UtCSID --pretrained_model ~/projects/def-stoeber/shrestha/Cellpose/models/sickle20230404  --diameter 0  --no_npy --save_txt --verbose --use_gpu
        done
done
