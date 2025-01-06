The generation of Fig. 3a requires three scripts in the 'UCYN-A_GitHub' folder, the 'UCYNA_Vary_Light.m', the 'UCYNA_Vary_Light_Size.m', and the 'Plot_Nfix.m'.

(1) First we need to run 'UCYNA_Vary_Light.m' to get the N2 fixation rate at different light intensities. This will also generate a light vs N2 fixation rate plot. It will take approximately 20 seconds to run. However, an increase in the step length of L_range from 0:1:800 to 0:20:800 will substantially reduce the runtime.

(2) To obtain the range of variations in N2 fixation rate at different light intensities, run 'UCYNA_Vary_Light_Size.m'. Make sure to keep the same step length of L_range as in (1).

(3) Then run 'Plot_Nfix.m' to generate the plot.

We also need 'viscocity_temperature_Jumars_1993.csv' file to run (1) and (2).
