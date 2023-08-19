# Frequency distribution with 100 bins and 40 morphological parameters
This section includes the calculation and plotting of frequency distribution with 100 bins for the 40 morphological parameters. 

The MATLAB program `FrequencyDistribution1D_100b40mp.m` is run in the cluster using the shell script `runFreq100b40mp.sh` using the following command (conservative estimate of time provided): 
```
sbatch --time=30:00:00 --array=1-4 runFreq30b40mp.sh
```
