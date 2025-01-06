# UCYN-A
Files related to the paper 'Nitrogen fixation mechanisms and distribution patterns of UCYN-A in the global ocean'


# Description of files
The folder 'UCYN-A_GitHub' contains different files.

(1) Matlab codes to generate Fig. 3a. The provided code will help to generate the variations in N2 fixation rate of UCYN-A at different light intensities.

(2) '.xlsx' files for the new experimental data used in Fig. 3 and Fig. 4.


# How to run the code
The generation of Fig. 3a requires three scripts in the 'UCYN-A_GitHub' folder, the 'UCYNA_Vary_Light.m', the 'UCYNA_Vary_Light_Size.m', and the 'Plot_Nfix.m'. 

(1) First we need to run 'UCYNA_Vary_Light.m' to get the N2 fixation rate at different light intensities. This will also generate a light vs N2 fixation rate plot. It will take approximately 20 seconds to run. However, an increase in the step length of L_range from 0:1:800 to 0:20:800 will substantially reduce the runtime.

(2) To obtain the range of variations in N2 fixation rate at different light intensities, run 'UCYNA_Vary_Light_Size.m'.

(3) Then run 'Plot_Nfix.m' to generate the plot. 

We also need 'viscocity_temperature_Jumars_1993.csv' file to run (1) and (2).


# Matlab download
To download Matlab, you will need to create a MathWorks account. Follow these steps:

(1) Open a MathWorks account at https://www.mathworks.com/login/.
(2) Once logged in, select "Get MATLAB" and then click on the "Start a Free Trial" button.
(3) You will be prompted to download Matlab, which is available for Windows, Linux, and Mac computers.
(4) Follow the on-screen instructions to install Matlab. The installation process may take up to 30 minutes.
(5) The free trial provides unlimited access for 30 days.


# Matlab version
The code was tested using Matlab R2022b.
