%Test Script to Run a PSO-based Matched-Filter search in the 
% positive chirp time space (S_P) or  negative chirp time space (S_N)

%Define paths
jsonlabpath = '/path/to/jsonlab-2.0';
ana_path = '/path/to/Accelerated-Network-Analysis/OptimizedPSOCodes';
psocodepath = '/path/to/Desktop/SDMBIGDAT19/CODES';

%Add paths
addpath(jsonlabpath);
addpath(ana_path);
addpath(psocodepath);


%Test 1
rungwpsosearch(1,'P');

%Test 2 
% rungwpsosearch(5,'N');

