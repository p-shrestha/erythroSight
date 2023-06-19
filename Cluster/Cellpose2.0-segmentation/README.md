# Segmentation of cells using Cellpose 2.0

The processed differential phase contrast (DPC) images were segmented using [Cellpose 2.0](https://www.nature.com/articles/s41592-022-01663-4), where a pre-trained neural network model (Cyto) was improved using a human-in-the-loop approach. The human-in-the-loop approach used 125 image segments from multiple donors with user annotations for each image segment and a recursive training approach. The annotations and training were fine-tuned to segment red blood cells (healthy/round and sickle shapes) and to avoid segmentation of cells at the edges of the images and overlapping cells.

## Nepal data segmentation
The different unique IDs or UIDs from the Nepal study are listed in `case_list` in ascending order of the size of the directories. The shell scripts `runCP-1.sh` to `runCP-4.sh` run Cellpose2.0 for different segments of the case_list at different processing times. The job arrays are submitted using `sbatch runCP-1.sh` to `sbatch runCP-4.sh` for processing the Nepal data. 

## Canada data segmentation
The different unique IDs or UIDs from the Canada study are listed in `case_listc` in ascending order of the size of the directories. The shell scripts `runCPc.sh` runs Cellpose2.0 for different segments case_listc. Note: Unlike the earlier implementation (Nepal data), the `--time` and `--array` information are not provided in the shell script, but rather entered when submitting the job, as follows:

```sbatch --time=16-10:00:00 --array=29 runCPc.sh
sbatch --time=14-10:00:00 --array=27-28 runCPc.sh
sbatch --time=12-12:00:00 --array=25-26 runCPc.sh
sbatch --time=12-00:00:00 --array=23-24 runCPc.sh
sbatch --time=11-00:00:00 --array=21-22 runCPc.sh
sbatch --time=9-20:00:00 --array=18-20 runCPc.sh
sbatch --time=9-05:00:00 --array=16-17 runCPc.sh
sbatch --time=8-01:00:00 --array=12-15 runCPc.sh
sbatch --time=6-14:00:00 --array=10-11 runCPc.sh
sbatch --time=6-09:00:00 --array=8-9 runCPc.sh
sbatch --time=6-00:00:00 --array=6-7 runCPc.sh
sbatch --time=4-15:00:00 --array=4-5 runCPc.sh
sbatch --time=4-00:00:00 --array=1-3 runCPc.sh

The time series data for Canada is processed using another script, which can be found in the TimeSeries folder. 
