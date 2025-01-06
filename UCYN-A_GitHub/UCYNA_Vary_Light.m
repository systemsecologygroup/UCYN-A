
clear all
clc
close all

% L = 300; % mumol/m^2/s
O2_E = 300; % mumol O2/L % O2
temp_current = 20; % Temperature
Fe = 1; % nmol Fe/L % Iron
P = 1; % mumol/L % Phosphate

%%%%%%%%%%%%% Independent parameters %%%%%%%%%%%%%

% max photosynthetic rate parameter
p_a = 5.5; 
p_b = 6.0;
p_c = 3.8; 

% light affinity parameter
a_L = 5;
b_L = 2; 

% Fe uptake rate parameter 
a_U_Fe = 6e-10; % nmol Fe/cell/d (Ward et al. 2012)
b_U_Fe = 0.6;  

% Fe affinity parameter
a_alpha_Fe = 2.5e-9; % L/cell/d  
b_alpha_Fe = 0.8; 

% P uptake rate parameter
a_U_P = 1.5e-10; % mumol P/cell/d
b_U_P = 0.77; 

% P affinity parameter
a_alpha_P = 3e-9; % L/cell/d 
b_alpha_P = 0.72;

% Max N2 fix rate
J_Nfix_max = 20; % 1*C_S; % 6; % fmol N/cell/d

% Cellular ratios
n_CN = 6.3; % mol C/mol N
n_O2C_photo = 1.4e-9; % mumol O2/fmol C % Freitas et al. 2020 (Photosynthetic quotient)
n_O2C_resp = 1e-9; % mumol O2/fmol C % Freitas et al. 2020 (Respiratory quotient)
n_FeC = 12*1e+3*1e-15; % nmol Fe/fmol C % Tuit et al., LO, 2004
n_PC = 1e-09*1/200; % mumol P/fmol C % Knapp et al. Aq Microbiol Ecol. 2012

% Respiration rates
R_0 = 0.1; % 0.04; % 1/d % Less than Chakraborty 2017 as host dont take DIN and UCYN-A dont do photosynthesis
R_L = 0.1; % 0.08; % fmol C/fmol C
R_N2 = 2.08; % 0.39; % fmol C/fmol N 
R_Fe = R_L/n_FeC; % 0.08; % fmol C/nmol Fe
R_P = R_L/n_PC; % 0.08; % fmol C/mumol P

% Diffusion coefficients 
D_O2_water = 2.12*1e-5; % McCabe and Laurent 1975  (1.38-2.12*10^(-5); % cm^2/s Broecker & Peng 1974)  (2*10^(-5) cm^2/s; Picioreanu et al 1997) %(diff coeff of O2 in water) 
eps_m_H = 1; % unitless (Inomura et al. 2017) %(diff. of O2 in cytoplasm relative to water)
eps_m_S = 1e-3; % 1.9e-3; % unitless (Inomura et al. 2017) %(diff. of O2 in cytoplasm relative to water)

% Temperature effect 
K = 273.15;

% For temperature dependent Q10 
a1 = 0.2;
a2 = 6;
a3 = 4;
a4 = 2.4;

% Q10 values
% Q10_J_L = 1.88; % Eppley 1972
Q10_J_L = a3 - a4*(1./(1+exp(-a1*temp_current))).^a2; % 2.2;%
Q10_R = 2; % Tait_2013 % 2; % Eppley 1972; Clarke 2003
Q10_U = 2.1; % (Ref in Maranon 2018) Raven JA, Geider RJ. Temperature and algal growth. New Phytol. 1988;110:441–61.

% Reference temperature 
ref_temp_DO2 = 25; % McCabe and Laurent 1975
ref_temp_J_L_max_H = 20; % assumption 
ref_temp_R_0 =20; % Mislan et al. 2014
ref_temp_R_L = 20; % Chakraborty et al. 2021
ref_temp_R_N2 = 28; % Grosskopf and LaRoche 2012
ref_temp_J_Fe_max = 20; % Ward 2012
ref_temp_J_P_max = 20; % ---

% Light attenuation rate
K_Z = 0.04; % 1/m % Letelier et al. 2004

% viscocity data        
vis_data = importdata('viscocity_temperature_Jumars_1993.csv');
vis_temp_dat = vis_data(:,1); % C
vis_dat = vis_data(:,2); % ﻿ g/cm/s

ref_vis_DO2 = interp1(vis_temp_dat,vis_dat,ref_temp_DO2,'linear','extrap'); 

% Range of investments 
theta_range = 0:0.001:1; % Fraction of fixed C transferred to UCYN-A
phi_range = 0:0.001:1; % Fraction of fixed N transferred to the host

rad_threshold = 0.005; % Threshold difference between the growth rate

L_range = 0:1:800; % Range of light intensity








%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% Vary Light %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%% Parameters dependent on size %%%%%%%%%%%%%%

% mass
diam_S = 0.6; % 0.65; % 0.8 % mum
rad_S = diam_S/2; % mum
vol_S = 4/3*pi*(rad_S)^3; % 0.27; % mum^3/cell

diam_H = 1.4; % 1.4; % 1.5; % 1.44; % mum
rad_H = diam_H/2; % mum
vol_H = 4/3*pi*(rad_H)^3; % 1.56; % mum^3

C_S = 10^(-0.363+0.863*log10(vol_S))*1e3/12; % fmol C/cell (% Verity et al. 1992)
C_H = 10^(-0.363+0.863*log10(vol_H))*1e3/12; % fmol C/cell (% Verity et al. 1992)

% max photosynthetic rate can be increases as the host does not take N
J_L_max_H = C_H*(p_a+log10(vol_H))/(p_b - p_c*log10(vol_H)+log10(vol_H)^2); % fmol C/cell/d 

% affinity for light
alpha_H = a_L*vol_H.^(2/3).*b_L*vol_H./(b_L*vol_H+a_L*vol_H.^(2/3));

% max Fe uptake rate
J_Fe_max_H = a_U_Fe*vol_H^b_U_Fe; % nmol Fe/cell/d
J_Fe_max_S = a_U_Fe*vol_S^b_U_Fe; % nmol Fe/cell/d

% affinity for Fe uptake
alpha_Fe_H = a_alpha_Fe*vol_H^b_alpha_Fe; % L/cell/d
alpha_Fe_S = a_alpha_Fe*vol_S^b_alpha_Fe; % L/cell/d

% max P uptake rate
J_P_max_H = a_U_P*vol_H^b_U_P; % mumol P/cell/d
J_P_max_S = a_U_P*vol_S^b_U_P; % mumol P/cell/d

% affinity for P uptake
alpha_P_H = a_alpha_P*vol_H^b_alpha_P; % L/cell/d
alpha_P_S = a_alpha_P*vol_S^b_alpha_P; % L/cell/d

%%%%%%%%%%%%%% Diffusion coefficients %%%%%%%%%%%%%%
Cw = 10*1e-7; % cm % (thickness of cell wall)
Lm_H = (1e-4*rad_H-Cw)/30; % cm % (thickness of plasma membrane)
Lm_S = (1e-4*rad_S-Cw)/30; % cm % (thickness of plasma membrane)
rc_H = 29*Lm_H; % cm %(thickness of cytoplasm)             
rc_S = 29*Lm_S; % cm %(thickness of cytoplasm)             

%%%%%%%%%%%%% Q10 %%%%%%%%%%%
J_L_max_H_Q10 = J_L_max_H*Q10_J_L^((temp_current-ref_temp_J_L_max_H)/10);
R_0_Q10 = R_0*Q10_R^((temp_current-ref_temp_R_0)/10);
R_L_Q10 = R_L*Q10_R^((temp_current-ref_temp_R_L)/10);
R_N2_Q10 = R_N2*Q10_R^((temp_current-ref_temp_R_N2)/10);
J_Fe_max_H_Q10 = J_Fe_max_H*Q10_U^((temp_current-ref_temp_J_Fe_max)/10);
J_Fe_max_S_Q10 = J_Fe_max_S*Q10_U^((temp_current-ref_temp_J_Fe_max)/10);
J_P_max_H_Q10 = J_P_max_H*Q10_U^((temp_current-ref_temp_J_P_max)/10);
J_P_max_S_Q10 = J_P_max_S*Q10_U^((temp_current-ref_temp_J_P_max)/10);
R_Fe_Q10 = R_Fe*Q10_R^((temp_current-ref_temp_R_L)/10); % same as R_L
R_P_Q10 = R_P*Q10_R^((temp_current-ref_temp_R_L)/10); % same as R_L

%%%%%%%%%%%%%% Diffusion temperature dependence  %%%%%%%%%%%%%%
vis = interp1(vis_temp_dat,vis_dat,temp_current,'linear','extrap'); 

D_O2_water_Q10 = D_O2_water*ref_vis_DO2*(K+temp_current)/(vis*(K+ref_temp_DO2));

K_O2_H = D_O2_water_Q10*1e-4*(rc_H+Lm_H)*eps_m_H/(rc_H*eps_m_H+Lm_H); % cm^2/s %(apparent diff of O2 in cell) 
D_H = 4*pi*1e-5*rad_H*1e-2*K_O2_H*24*60*60; % dm^3/d = L/d/cell %(oxygen flux in cell)
    
K_O2_S = D_O2_water_Q10*1e-4*(rc_S+Lm_S)*eps_m_S/(rc_S*eps_m_S+Lm_S); % cm^2/s %(apparent diff of O2 in cell) 
D_S = 4*pi*1e-5*rad_S*1e-2*K_O2_S*24*60*60; % dm^3/d = L/d/cell %(oxygen flux in cell)

%%%%%%%%%%%%%% Memory allocation %%%%%%%%%%%%%%

Nfix = zeros(1,length(L_range));

for i1 = 1:length(L_range)
    
    L = L_range(i1)

    %%%%%%%%%%%%%% photosynthesis and Fe uptake rate  %%%%%%%%%%%%%%
    
    J_L = J_L_max_H_Q10*(1-exp(-alpha_H*L/J_L_max_H_Q10)); % exp(-beta*L/J_L_max); % fmol C/cell-H/day     
    
    J_Fe_H_all = J_Fe_max_H_Q10*alpha_Fe_H*Fe/(J_Fe_max_H_Q10+alpha_Fe_H*Fe); % nmol Fe/cell/d
    J_Fe_H = 0.85*J_Fe_H_all; 
    J_Fe_S = 0.15*J_Fe_H_all; % J_Fe_max_S_Q10*alpha_Fe_S*Fe/(J_Fe_max_S_Q10+alpha_Fe_S*Fe); % nmol Fe/cell/d
    
    J_P_H_all = J_P_max_H_Q10*alpha_P_H*P/(J_P_max_H_Q10+alpha_P_H*P); % mumol P/cell/d
    J_P_H = 0.85*J_P_H_all;
    J_P_S = 0.15*J_P_H_all; % J_P_max_S_Q10*alpha_P_S*P/(J_P_max_S_Q10+alpha_P_S*P); % mumol P/cell/d

    %%%%%%%%%%%%%% Growth rate calculation %%%%%%%%%%%%%%
    
    % ------- Host -------- %
                             
    J_HC_bar_temp = max(0,J_L - R_0_Q10*C_H - R_L_Q10*J_L - R_Fe_Q10*J_Fe_H_all - R_P_Q10*J_P_H_all); % fmol C/cell/day
    Resp_H = R_0_Q10*C_H + R_L_Q10*J_L + R_Fe_Q10*J_Fe_H_all + R_P_Q10*J_P_H_all; % fmol C/cell/day
    
    O2_H = max(0, D_H*O2_E + n_O2C_photo*J_L - n_O2C_resp*Resp_H)/(D_H+D_S);
    
    div_rate_H_temp = zeros(length(theta_range),length(phi_range));

    for i = 1:length(theta_range)
        theta = theta_range(i);
        
        for j = 1:length(phi_range)
            phi = phi_range(j);
            
            % ------- Symbiont -------- %
                
            J_SN_temp = phi*J_Nfix_max; % fmol N/cell/day
    
            R1_temp = R_0_Q10*C_S + R_N2_Q10*J_SN_temp; % fmol C/cell/day    
                
            R_O2_temp = max(0,D_S*O2_H/n_O2C_resp - R1_temp ); % fmol C/cell/day % remaining amount of o2 to remove

            if R_O2_temp > 0
                phi = 0;
                J_SN_temp = 0; % fmol N/cell/day
                R1_temp = R_0_Q10*C_S; % fmol C/cell/day    
            end

            Resp_S_temp = R1_temp;% + R_O2_temp;
    
            J_SC_temp = max(0, theta*J_HC_bar_temp - Resp_S_temp); % fmol C/cell/day
            J_SG_temp = min(min(min(J_SC_temp, n_CN*J_SN_temp), J_Fe_S/n_FeC), J_P_S/n_PC); % fmol C/cell/day
            div_S_temp = max(0,J_SG_temp/C_S); % 1/day
            rho_N_temp = max(0,J_SN_temp - J_SG_temp/n_CN); % fmol N/cell/d
                    
            % ------- Host -------- %
                             
            J_HC_temp = (1-theta)*J_HC_bar_temp; % fmol C/cell/day
            J_HN_temp = rho_N_temp; % fmol N/cell/day
            J_HG_temp = min(min(min(J_HC_temp, n_CN*J_HN_temp), J_Fe_H/n_FeC), J_P_H/n_PC); % fmol C/cell/day
            div_H_temp = J_HG_temp/C_H; % 1/day
            
            % ------- Measure -------- %
                
            if (div_S_temp > 1e-3) && (div_H_temp > 1e-3) 
                if abs(div_S_temp - div_H_temp) < rad_threshold
                    div_rate_H_temp(i,j) = div_H_temp;
                end
            else
                div_rate_H_temp(i,j) = NaN;
            end
        
        end
    end
    
    [val,ind] = max(div_rate_H_temp(:));
    [ind_row, ind_col] = ind2sub(size(div_rate_H_temp),ind);
    
    Nfix(i1) = phi_range(ind_col)*J_Nfix_max;
    
    ind = [];
    ind_row = [];
    ind_col = [];
    theta_val = [];
end

save L_range.mat L_range
save Nfix.mat Nfix


%% Plot N2 fixation rate 

plot(L_range,Nfix,'Color','b','LineWidth',2);
hold on;

%% Set plot requirements

xlim([0 800])
xlabel('PAR (\mumol m^{-2} s^{-1})')
xticks(0:200:800)
xticklabels({'0','200','400','600','800'})

ylim([0 10])
ylabel('N_2 fixation rate (fmol N cell^{-1} d^{-1})')
yticks(0:2:10)
yticklabels({'0','2','4','6','8','10'})

set(gca,'TickLength',[0.015, 0.01])
set(gca,'FontSize',13)

box on
hold off;

