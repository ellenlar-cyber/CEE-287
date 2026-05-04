%% HW 3 Part 7
clc; clear;close all;

% BRBF parameters from ASCE 7-22 Table 12.2-1
R  = 8;
Cd = 5;
Ie = 1.0;   % Risk Category II
T  = 0.6;   % fundamental period (s)
g  = 981;   % cm/s^2

%% Load ground motions
data     = readmatrix('Ground_motions_assignment_1.csv');
t_record = data(:,1);
dt       = t_record(2) - t_record(1);
g        = 981;

% Columns: 1=time, 2=Parking, 3=SLAC-1, 4=SLAC-2, 5=VA-1, 6=VA-2
ag_SLAC1 = data(:,3);
ag_SLAC2 = data(:,4);
ag_VA1   = data(:,5);
ag_VA2   = data(:,6);


%% Final scale factors DBE
SF_SLAC_DBE = 2.8634;
SF_VA_DBE   = 2.7373;


%% Final scale factors MCE
SF_SLAC_MCE = 4.2972;
SF_VA_MCE   = 4.1080;

%% Scale ground motion records for DBE
ag_SLAC1_sc_d = ag_SLAC1 * SF_SLAC_DBE;
ag_SLAC2_sc_d = ag_SLAC2 * SF_SLAC_DBE;
ag_VA1_sc_d   = ag_VA1   * SF_VA_DBE;
ag_VA2_sc_d   = ag_VA2   * SF_VA_DBE;

%% Scale ground motion records for MCE
ag_SLAC1_sc_m = ag_SLAC1 * SF_SLAC_MCE;
ag_SLAC2_sc_m = ag_SLAC2 * SF_SLAC_MCE;
ag_VA1_sc_m  = ag_VA1   * SF_VA_MCE;
ag_VA2_sc_m   = ag_VA2   * SF_VA_MCE;

%% Compute individual response spectra of scaled records (2% damping)
z  = 0.02;
wn = 2*pi/T;

% Calculate Cy 
% Cy = Sa design/ R
Sa_design_dbe = 1.4800;
Sa_design_mce = 2.2200; 

Cy_dbe = Sa_design_dbe/R;
Cy_mce = Sa_design_mce/R;

%% DBE Values using inputs and non linear response function developed in previous hw
    [~, ~, ~, ~, Sd_in_d_S1, ~] = SDOF_Response_NL_1...
        (T, z, ag_SLAC1_sc_d, dt, 0, 0, Cy_dbe, 'Cy', 'linear');

    [~, ~, ~, ~, Sd_in_d_S2, ~] = SDOF_Response_NL_1...
        (T, z, ag_SLAC2_sc_d, dt, 0, 0, Cy_dbe, 'Cy', 'linear');

    [~, ~, ~, ~, Sd_in_d_V1, ~] = SDOF_Response_NL_1...
        (T, z, ag_VA1_sc_d, dt, 0, 0, Cy_dbe, 'Cy', 'linear');

    [~, ~, ~, ~, Sd_in_d_V2, ~] = SDOF_Response_NL_1...
        (T, z, ag_VA2_sc_d, dt, 0, 0, Cy_dbe, 'Cy', 'linear');

    Sd_d_mean = (Sd_in_d_S1 + Sd_in_d_S2 + Sd_in_d_V1 + Sd_in_d_V2) / 4;

    Sd_d_max = max( [Sd_in_d_S1, Sd_in_d_S2, Sd_in_d_V1, Sd_in_d_V2 ]);

    
%% MCE values using inputs and non linear response function developed in previous hw
    [~, ~, ~, ~, Sd_in_m_S1, ~] = SDOF_Response_NL_1...
        (T, z, ag_SLAC1_sc_m, dt, 0, 0, Cy_mce, 'Cy', 'linear');

    [~, ~, ~, ~, Sd_in_m_S2, ~] = SDOF_Response_NL_1...
        (T, z, ag_SLAC2_sc_m, dt, 0, 0, Cy_mce, 'Cy', 'linear');

    [~, ~, ~, ~, Sd_in_m_V1, ~] = SDOF_Response_NL_1...
        (T, z, ag_VA1_sc_m, dt, 0, 0, Cy_mce, 'Cy', 'linear');

    [~, ~, ~, ~, Sd_in_m_V2, ~] = SDOF_Response_NL_1...
        (T, z, ag_VA2_sc_m, dt, 0, 0, Cy_mce, 'Cy', 'linear');

    Sd_m_mean = (Sd_in_m_S1 + Sd_in_m_S2 + Sd_in_m_V1 + Sd_in_m_V2) / 4;

    Sd_m_max = max( [Sd_in_m_S1, Sd_in_m_S2, Sd_in_m_V1, Sd_in_m_V2 ]);

fprintf('\nSd: DBE Spectrum, 2%% Damping\n');
fprintf('Sd_DBE_S1 = %.4f cm\n', Sd_in_d_S1);
fprintf('Sd_DBE_S2 = %.4f cm\n', Sd_in_d_S2);
fprintf('Sd_DBE_VA1 = %.4f cm\n', Sd_in_d_V1);
fprintf('Sd_DBE_VA2 = %.4f cm\n', Sd_in_d_V2);
fprintf('Sd_DBE_mean = %.4f cm\n', Sd_d_mean);
fprintf('Sd_DBE_max = %.4f cm\n', Sd_d_max);

fprintf('\nSd: MCE Spectrum, 2%% Damping\n');
fprintf('Sd_MCE_S1 = %.4f cm\n', Sd_in_m_S1);
fprintf('Sd_MCE_S2 = %.4f cm\n', Sd_in_m_S2);
fprintf('Sd_MCE_VA1 = %.4f cm\n', Sd_in_m_V1);
fprintf('Sd_MCE_VA2 = %.4f cm\n', Sd_in_m_V2);
fprintf('Sd_MCE_mean = %.4f cm\n', Sd_m_mean);
fprintf('Sd_MCE_max = %.4f cm\n', Sd_m_max);

    