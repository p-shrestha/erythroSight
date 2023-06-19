# Segmentation for time series data for Canada

The different unique IDs or UIDs from the time series data of the Canada study are listed in `caseListTimeSeries` in ascending order of the size of the directories. The shell scripts `runCPTimeSeries.sh` runs Cellpose2.0 for different segments `caseListTimeSeries`. The jobs can be submitted using the following commands: 

```
sbatch --time=12:00:00 --array=1 runDPCTimeSeries.sh 
sbatch --time=20:00:00 --array=2-9 runDPCTimeSeries.sh 
sbatch --time=1-15:00:00 --array=10-12 runDPCTimeSeries.sh 
sbatch --time=2-11:00:00 --array=13-14 runDPCTimeSeries.sh 
sbatch --time=3-00:00:00 --array=15 runDPCTimeSeries.sh
```
