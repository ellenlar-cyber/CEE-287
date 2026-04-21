%% HW 3 Part 5 - Compute Displacement Demands

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

%% DBE
    [~,~,~,Sd_S1_d,~,~,~,~] = SDOF_Response(T, z, ag_SLAC1_sc_d, dt, 0, 0);


    [~,~,~,Sd_S2_d,~,~,~,~] = SDOF_Response(T, z, ag_SLAC2_sc_d, dt, 0, 0);


    [~,~,~,Sd_V1_d,~,~,~,~] = SDOF_Response(T, z, ag_VA1_sc_d, dt, 0, 0);
 

    [~,~,~,Sd_V2_d,~,~,~,~] = SDOF_Response(T, z, ag_VA2_sc_d, dt, 0, 0);
  

    Sd_mean_DBE = (Sd_S1_d + Sd_S2_d + Sd_V1_d + Sd_V2_d) / 4;

    Sd_max_DBE = max([Sd_S1_d, Sd_S2_d, Sd_V1_d, Sd_V2_d]);

%% MCE
    [~,~,~,Sd_S1_m,~,~,~,~] = SDOF_Response(T, z, ag_SLAC1_sc_m, dt, 0, 0);

    [~,~,~,Sd_S2_m,~,~,~,~] = SDOF_Response(T, z, ag_SLAC2_sc_m, dt, 0, 0);

    [~,~,~,Sd_V1_m,~,~,~,~] = SDOF_Response(T, z, ag_VA1_sc_m, dt, 0, 0);

    [~,~,~,Sd_V2_m,~,~,~,~] = SDOF_Response(T, z, ag_VA2_sc_m, dt, 0, 0);

    Sd_mean_MCE = (Sd_S1_m + Sd_S2_m + Sd_V1_m + Sd_V2_m) / 4;

    Sd_max_MCE = max([Sd_S1_m, Sd_S2_m, Sd_V1_m, Sd_V2_m]);

% Eq. 12.8-16: DBE displacement
delta_DBE_S1 = ((Cd / R) / Ie) * Sd_S1_d;
delta_DBE_S2  = ((Cd / R) / Ie) * Sd_S2_d;
delta_DBE_V1 = ((Cd / R) / Ie) * Sd_V1_d;
delta_DBE_V2 = ((Cd / R) / Ie) * Sd_V2_d;
delta_DBE_mean = ((Cd / R) / Ie) * Sd_mean_DBE;
delta_DBE_max = ((Cd / R) / Ie) * Sd_max_DBE;


% Eq. 12.8-17: MCE displacement ( DO NOT SCALE BY 1.5 BECAUSE ALREADY USING
% MCE SPECTRUM, NOT CONVERTING FROM DBE )
delta_MCE_S1 = 1 / Ie * Sd_S1_m;
delta_MCE_S2 = 1  / Ie * Sd_S2_m;
delta_MCE_V1 = 1 / Ie * Sd_V1_m;
delta_MCE_V2 = 1 / Ie * Sd_V2_m;
delta_MCE_mean = 1 / Ie * Sd_mean_MCE;
delta_MCE_max = 1 / Ie * Sd_max_MCE;

fprintf('delta_DBE_S1 = %.4f cm\n', delta_DBE_S1);
fprintf('delta_DBE_S2 = %.4f cm\n', delta_DBE_S2);
fprintf('delta_DBE_VA1 = %.4f cm\n', delta_DBE_V1);
fprintf('delta_DBE_VA2 = %.4f cm\n', delta_DBE_V2);
fprintf('delta_DBE_mean = %.4f cm\n', delta_DBE_mean);
fprintf('delta_DBE_max = %.4f cm\n', delta_DBE_max);


fprintf('delta_MCE_S1 = %.4f cm\n', delta_MCE_S1);
fprintf('delta_MCE_S2 = %.4f cm\n', delta_MCE_S2);
fprintf('delta_MCE_VA1 = %.4f cm\n', delta_MCE_V1);
fprintf('delta_MCE_VA2 = %.4f cm\n', delta_MCE_V2);
fprintf('delta_MCE_mean = %.4f cm\n', delta_MCE_mean);
fprintf('delta_MCE_max = %.4f cm\n', delta_MCE_max);
