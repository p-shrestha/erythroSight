# 1D frequency characterization for 30 bins and 16 morphological parameters 

This section includes programs, plots, data and results for 30 bins and 16 morphological parameters (Area,Perimeter,Angle,Major,Minor,Height,Width,Feret,MinFeret,Circularity,Roundness,Solidity,Eccentricity,ESF,ElongationFmF,Convexity).

All the `.sh` shell scripts run the MATLAB program `FrequencyDistribution1D_30b16mp.m` 

The jobs can be submitted using the following commands in the cluster: 
```
sbatch runFreq1D-AA.sh
sbatch runFreq1D-ABeta.sh
sbatch runFreq1D-AS.sh
sbatch runFreq1D-SBeta.sh
sbatch runFreq1D-SS.sh
sbatch runFreq1D-N.sh
sbatch runFreq1D-SCT.sh
sbatch runFreq1D-SCD.sh
```
