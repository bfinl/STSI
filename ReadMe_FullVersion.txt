This README document contains the code and de-identified human testing data for the Spatial-Temporal-Spectral Source Imaging (STSI) algorithm.

The source code and testing data are provided as a service to the scientific community and may be used for any non-commercial purposes under a CC-BY-NC-SA-4.0 license. Users should use the provided codes or data at their own risk. A copy of the CC-BY-NC-SA-4.0 license is provided along with this program. If not, see https://creativecommons.org/licenses/by-nc-sa/4.0/legalcode.

As part of the license, please cite the paper below if you use the algorithms or results in your research:

Jiang X, Cai Z, Gonsisko C, Worrell G, He B. Mapping epileptogenic brain using a unified spatial–temporal–spectral source imaging framework. Proc Natl Acad Sci U S A. 2025;122:e2510015122. doi:10.1073/pnas.2510015122

This work was funded by the National Institutes of Health grants NS096761, NS127849, NS131069, NS124564, and EB029365 (B.H.).

•	Software requirements
The code provided is based on MATLAB, and mainly tested in Windows platform with MATLAB version 2024b. If problems occur when running the code, it is likely a compatibility issue, e.g. the operating system or the MATLAB version. 

Toolboxes
Please ensure the following toolboxes are available in your MATLAB environment:
•	MATLAB Wavelet Toolbox: available through the MATLAB installation (licensed).
•	Tensorlab: obtain from the official website: https://tensorlab.net/.
•	EEGLAB: obtain from the official website: https://eeglab.org/download/.
•	FieldTrip: obtain from the official website: https://www.fieldtriptoolbox.org/download/.
To avoid potential issues, one can remove the path of other versions of these toolboxes from the MATLAB path. Please put Tensorlab, EEGLAB, and FieldTrip toolbox into the 'ToolBox_External' folder for the code to load.

•	Setting up the path
The startup_patient.m is the script to set up relevant paths. The current version assumes the working directory is E:\STSI_CodesAndNotes\STSI_Codes\, if your drive is different, please edit the startup script: ...\Codes\startup_patient.m
And change the variable MainPath = 'E:\STSI_CodesAndNotes\STSI_Codes\' to your own directory.

•	Explanation of the subfolders
Codes: the example codes from one patient. The functions and code related to the algorithms are stored in Folder ‘RelevantFunctions’. 
- Area.m: the function to calculate triangle area.
- CortexViewInit.m: initialize the cortex viewing angle.
- norms.m: calculates the norm of a vector (or a set of vectors).
- cvx_check_dimension.m: Verifies that the input is valid dimension, used as part of the function called ‘norms’.
- cvx_default_dimension.m: Default dimension for SUM, MAX, etc, used as part of the function called ‘norms’.
- Find_Patch.m: find continuous source imaging results patches for display
- FISTA_ADMM_IRES_Tensor.m: the projection to hyperellipsoid method during iterative solving of STSI.
- trisurf_customized.m: display triangulated surface such as the cortex.
Figures: stores source imaging results figures
Grid_Location_Parameters: stores the head model, including leadfield matrices, neighboring relationships, electrode positions, etc.
Seizures: stores the seizure data, the data format is a N_channel by N_sample matrix.
Spikes: stores the HFO riding spike data, the data format is a N_channel by N_sample matrix.

•	How to run the code
HFO source imaging, navigate to:
STSI_CodesAndNotes\STSI_Codes\Codes\ Example_Patient1_Event1_pHFO.m. 
Spike source imaging:
STSI_CodesAndNotes\STSI_Codes\Codes\ Example_Patient1_Event1_pSpike.m
Seizure source imaging: 
STSI_CodesAndNotes\STSI_Codes\Codes\ Example_Patient1_Seizure1.m

Each code should be click-and run, the results will be Figures showing the source imaging results.

Due to the sensitive nature of epilepsy patient data, the patient data is de-identified, and the exact surgical resection site for the example patient cannot be shared. Instead, we indicate the qualitative location as the left temporal cortex. This patient achieved an ILAE Class I outcome.
