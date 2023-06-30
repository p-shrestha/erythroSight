# 1D Characterization of time series data 

Programs to plot frequency distribution and empirical cumulative distribution for time series data from Canada.

The MATLAB program `FrequencyDistribution1DTimeSeries.m` reads morphological parameter data and creates and saves the plots. The `runFreqTimeSeries.sh` shell script runs the job-array in the cluster, and can be submitted using command line code as follows: 
`sbatch --time=02:00:00 --array=1-15 runFreqTimeSeriesTiled.sh`
