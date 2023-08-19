# Frequency distribution with 30 bins and 40 morphological parameters
This section includes the calculation and plotting of frequency distribution with 30 bins for the 40 morphological parameters. 

This implementation of the MATLAB program `FrequencyDistribution1D_30b40mp.m` includes sections that run the program `normalizeBasicParam.m` in the beginning. If the basic parameters are already created and saved, the first three sections can be commented, and the next section uncommented. 

The MATLAB program `FrequencyDistribution1D_30b40mp.m` is run in the cluster using the shell script `runFreq30b40mp.sh` using the following command (conservative estimate of time provided): 
```
sbatch --time=30:00:00 --array=1-4 runFreq30b40mp.sh
```
