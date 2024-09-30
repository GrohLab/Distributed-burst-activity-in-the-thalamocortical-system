#### This folder contains all necessary Matlab files, used for obtaining the relevant figure data. The data is accessed through the respective user directory path to the downloaded data files.

**Before executing the scripts, include the folder path in the MATLAB environment by using the command: addpath(genpath('.\Scripts')).**

#### Necessary MATLAB Toolboxes (to be installed before executing the scripts)
	I) 'Signal Processing Toolbox'
	II) 'Deep Learning Toolbox'
	III) 'Statistics and Machine Learning Toolbox'
	IV) 'Communications Toolbox'
	V) 'Image Processing Toolbox'

| File | Description |
|--------------------------------------------|------------------------------------------------------------------------------------------------------------------------------|
| userDataPath.m		             | Sets the path to the data (assuming the default Documents dir as a default). Run before running other scripts.		    |
| singleUnitBurstinessAnalysis.m             | Analyzes the burstiness and firing rates of single units, generating plots showing PSTHs (Fig. 2b,e).   			    |
| ApertureResponseTypes.m                    | Analyzes burst indices from different brain areas and conditions (Fig. 3b).						    |
| ApertureResponseTypes_Comparison.m         | Compares the burst biases of individual cells from different brain areas (Fig. 3c, 5a).		                   	    |
| ApertureResponseTypes_StageProgression.m   | Visualizes the progression of burst bias and the proportion of burst-responding cells as a function of learning (Fig. 5b,c). |
| DecodingWithIncreasingUnitNum_Plotting.m   | Visualizes the performance of decoding analysis with increasing number of units for different recording areas (Fig. 6a-d).   |
