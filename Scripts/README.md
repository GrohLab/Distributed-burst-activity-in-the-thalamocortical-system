#### This folder contains all necessary Matlab files, used for obtaining the relevant figure data. The data is accessed through the respective user directory path to the downloaded data files.

**Before executing the scripts, include the folder path in the MATLAB environment by using the command: addpath(genpath('.\Scripts')).**

#### Necessary MATLAB Toolboxes (to be installed before executing the scripts)
	I) 'Signal Processing Toolbox'
	II) 'Deep Learning Toolbox'
	III) 'Statistics and Machine Learning Toolbox'
	IV) 'Communications Toolbox'
	V) 'Image Processing Toolbox'

| File | Description |
|----------|----------|
| StimulusResponse_Analysis.m                | Get the StimulusResponse variable                                                     |
| DE_Salience_function.m                     | Get the waveforms_all.mat and allResponses.mat variables                              |
| getBurstiness.m                            | Get the BurstinessData.mat variable                                                   |
| burstsUponTrigger.m                        | Get the ResponsePattern variable							     |
| ApertureResponseTypes.m                    | Get the responseTypes variable and generate the burstiness scatters (Fig. 3b)	     |	
| ApertureResponseTypes_Comparison.m         | For comparing the burstiness bias of individual cells (Fig. 3c, 5a)		     |
| ApertureResponseTypes_StageProgression.m   | To calculate the burstiness stage progression (Fig. 5b,c)			     |
| DecodingWithIncreasingUnitNum_Analysis.m   | To calculate the decoding accuracy (Fig. 6a-d)					     |
| DecodingWithIncreasingUnitNum_Plotting.m   | To plot the decoding accuracy (Fig. 6a-d)					     |
| SingleCellApertureReactivity.m             | Necessary for running the singleUnitBurstinessAnalysis.m script			     |
| singleUnitBurstinessAnalysis.m             | For population PSTHs (Fig. 2b,e)							     |
