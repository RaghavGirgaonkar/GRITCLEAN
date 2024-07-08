This file contains instructions to run codes for a PSO-based matched filtering search and to generate
relevant plots.

Any queries can be directed to raghav.girgaonkar@gmail.com
##########################################################################
1. Dependencies

i. These codes require the JSONLAB package for MATLAB (https://www.mathworks.com/matlabcentral/fileexchange/33381-jsonlab-a-toolbox-to-encode-decode-json-files)
ii. The Parallel Computing Toolbox and Signal Processing Toolbox for MATLAB also need to be installed. 
iii. Clone the repositories Accelerated-Network-Analysis (https://github.com/RaghavGirgaonkar/Accelerated-Network-Analysis)
    and SDMBIGDAT19 (https://github.com/mohanty-sd/SDMBIGDAT19). 
iv. Add paths to the above three packages, specifically to 'Accelerated-Network-Analysis/OptimizedPSOCodes' and 
     'SDMBIGDAT19/CODES' to MATLAB and in "test_rungwpsosearch.m"

###########################################################################
2. Downloading HDF5 files

All HDF5 files can be downloaded using the GWOSC website (https://gwosc.org)
The HDF5 files for the dataset can be downloaded in 2 ways.

i. The URLs for all files are present in "urls.txt". Simply run "wget -i urls.txt" in the terminal command line to download all.
ii. The relevant information of all files in present in "hdf5files.json" with a url provided for each file. 
    One can download the file by copy-pasting this link in a web browser.

iii. Additionally, the json files directly downloaded from GWOSC are provided along with this dataset.
     These are "O1file.json", "O2files.json" etc. One could also download the corresponding files using the urls
     given in these json files.

A sample file "H-H1_GWOSC_O3a_4KHZ_R1-1243394048-4096.hdf5" is provided with this dataset. 


###########################################################################
3. Running a PSO-based matched-filtering search

i. A PSO-based search can be launched using the "rungwpsosearch.m" script.
   The inputs for this script are the segment number and the type of search 
   "P" for an S_P search, "N" for an S_N search. 

ii. The information about segments are provided in the "segments.json" file.
    This file contains the GPS start and end times, the HDF5 file used as well as the training indices 
    used to estimate the PSD.

iii. The file "GVSsegPSDtrainidxs.mat" contains all training indices required for estimating the Power Spectral Density for 
     each segment. This file should be in the same directory as the codes. 

iv. The parameters need for running PSO, such as the number of independent runs, number of iterations along with
    the parameters such as the sampling frequency, the search range in the chirp-time space can be changed if need via the
    json files ("signal.json", "pso.json" etc). For more details on how to change these parameters, the user can refer to the
    "usermanual.pdf" file in the Accelerated-Network-Analysis GitHub repository. 

A test script "test_rungwpsosearch.m" has been provided. This can run PSO-based searches on the provided hdf5 file 
"H-H1_GWOSC_O3a_4KHZ_R1-1243394048-4096.hdf5". The segment numbers that can be provided for this file to the "rungwpsosearch.m"
script are from 1 to 9. To run the search on other segments, make sure to have the correct files downloaded. 

To run this script simply run 'test_rungwpsosearch' in MATLAB.  

############################################################################# 
4. Dataset and Generating plots

i. The directory 'Output/' contains text files having estimated parameters for all candidates across all segments.
   These files also include the estimates for different signal injections (real, low and high SNR).
   Each file, bar the noise files, has 4 columns, the columns provide the parameters as follows:
   Segment Number   Tau_0   Tau_1.5      TOA    SNR
    
   
ii. The script "plot_candidates.m" can generate all major scatter plots as given in the paper using the dataset of all 
    candidates provided. Simply run this script in MATLAB. 

iii. The plotting script has a few changeable parameters, namely, the detection threshold, the complex mass threshold
     and the Delta TOA and Delta SNR thresholds. To reproduce the plots in the paper, the following thresholds should be used
     
     a.) detection_threshold = 9, complexmass_threshold = 0.9, deltatoa_threshold = 0.15 and deltasnr_threshold = 0.1
     b.) detection_threshold = 8, complexmass_threshold = 0.9, deltatoa_threshold = 0.15 and deltasnr_threshold = 0.07


iv. To reproduce numbers from Table II of the paper, the user can change the complexmass_threshold values as given in the table.
    The plotting script displays the relevant statistics.

##############################################################################
All codes have developed and tested on MATLAB R2022 and MATLAB R2023.
############################################################################## 
