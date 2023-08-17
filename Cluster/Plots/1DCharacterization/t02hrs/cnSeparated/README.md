# Processing t02hrs data for Canada and Nepal data separated

The data is processed according to the following groups: 
1. Canada
  - AA
  - ABeta
  - AS
  - SBeta
  - SS  
2. Nepal
  - N
  - SCT
  - SCD

## FrequencyDistribution1Dp.m
This program primarily creates and saves `morphParamData.mat`, `generalInfo.mat` and `histData.mat` for each group. Importantly, the `morphParamData.mat` which contains the `morphParam` variables including the morphological parameters of each cell at t02hrs is saved in MATLAB readable format here (read from .csv files). The frequency distribution saved in `histData` is not critical and is calculated elsewhere for different configurations (e.g. 30b16mp) and combinations (cnCombined). 

The shell scripts are run using the following commands: 
```
sbatch Freq1D-AA.sh
sbatch Freq1D-ABeta.sh
sbatch Freq1D-AS.sh
sbatch Freq1D-SBeta.sh
sbatch Freq1D-SS.sh
sbatch Freq1D-N.sh
sbatch Freq1D-SCT.sh
sbatch Freq1D-SCD.sh
```
