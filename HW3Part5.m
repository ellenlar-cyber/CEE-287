%% HW 3 Part 5 - Compute Displacement Demands

clc; clear;close all;
%% Questions: should the delta e be different, should it be computed based on the 5% displacements? do they use the same Sa? should this be based on Sa or Sd?
%% Do we need a plot? or is just for T struct?
%% do we use the scales ground motions?


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

%% UPDATE IF NEEDED!!!!!!
%% Final scale factors DBE
SF_SLAC_DBE = 2.8634;
SF_VA_DBE   = 2.7373;

%% UPDATE IF NEEDED!!!!!!
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

%% DBE
    [~,~,~,Sd,~,~,~,~] = SDOF_Response(T, z, ag_SLAC1_sc_d, dt, 0, 0);
    Sa_S1_d = wn^2 * Sd / g;

    [~,~,~,Sd,~,~,~,~] = SDOF_Response(T, z, ag_SLAC2_sc_d, dt, 0, 0);
    Sa_S2_d = wn^2 * Sd / g;

    [~,~,~,Sd,~,~,~,~] = SDOF_Response(T, z, ag_VA1_sc_d, dt, 0, 0);
    Sa_V1_d = wn^2 * Sd / g;

    [~,~,~,Sd,~,~,~,~] = SDOF_Response(T, z, ag_VA2_sc_d, dt, 0, 0);
    Sa_V2_d = wn^2 * Sd / g;

    Sa_mean_DBE = (Sa_S1_d + Sa_S2_d + Sa_V1_d + Sa_V2_d) / 4;

    Sa_max_DBE = max([Sa_S1_d, Sa_S2_d, Sa_V1_d, Sa_V2_d]);

%% MCE
    [~,~,~,Sd,~,~,~,~] = SDOF_Response(T, z, ag_SLAC1_sc_m, dt, 0, 0);
    Sa_S1_m = wn^2 * Sd / g;

    [~,~,~,Sd,~,~,~,~] = SDOF_Response(T, z, ag_SLAC2_sc_m, dt, 0, 0);
    Sa_S2_m = wn^2 * Sd / g;

    [~,~,~,Sd,~,~,~,~] = SDOF_Response(T, z, ag_VA1_sc_m, dt, 0, 0);
    Sa_V1_m = wn^2 * Sd / g;

    [~,~,~,Sd,~,~,~,~] = SDOF_Response(T, z, ag_VA2_sc_m, dt, 0, 0);
    Sa_V2_m = wn^2 * Sd / g;

    Sa_mean_MCE = (Sa_S1_m + Sa_S2_m + Sa_V1_m + Sa_V2_m) / 4;

    Sa_max_MCE = max([Sa_S1_m, Sa_S2_m, Sa_V1_m, Sa_V2_m]);


% Calculate delta e based on Sa values
% delta_e = (Sa/(R/Ie)) * g*T^2/(4*pi^2)
delta_e_d_S1 = (Sa_S1_d / (R/Ie)) * g * T^2 / (4*pi^2);
delta_e_d_S2 = (Sa_S2_d / (R/Ie)) * g * T^2 / (4*pi^2);
delta_e_d_V1 = (Sa_V1_d / (R/Ie)) * g * T^2 / (4*pi^2);
delta_e_d_V2 = (Sa_V2_d / (R/Ie)) * g * T^2 / (4*pi^2);
delta_e_d_mean = (Sa_mean_DBE / (R/Ie)) * g * T^2 / (4*pi^2);
delta_e_d_max = (Sa_max_DBE / (R/Ie)) * g * T^2 / (4*pi^2);


% Calculate delta e based on Sa values
% delta_e = (Sa/(R/Ie)) * g*T^2/(4*pi^2)
delta_e_m_S1 = (Sa_S1_m / (R/Ie)) * g * T^2 / (4*pi^2);
delta_e_m_S2 = (Sa_S2_m / (R/Ie)) * g * T^2 / (4*pi^2);
delta_e_m_V1 = (Sa_V1_m / (R/Ie)) * g * T^2 / (4*pi^2);
delta_e_m_V2 = (Sa_V2_m / (R/Ie)) * g * T^2 / (4*pi^2);
delta_e_m_mean = (Sa_mean_MCE / (R/Ie)) * g * T^2 / (4*pi^2);
delta_e_m_max = (Sa_max_MCE / (R/Ie)) * g * T^2 / (4*pi^2);


% Eq. 12.8-16: DBE displacement
delta_DBE_S1 = (Cd / Ie) * delta_e_d_S1;
delta_DBE_S2  = (Cd / Ie) * delta_e_d_S2;
delta_DBE_V1 = (Cd / Ie) * delta_e_d_V1;
delta_DBE_V2 = (Cd / Ie) * delta_e_d_V2;
delta_DBE_mean = (Cd / Ie) * delta_e_d_mean;
delta_DBE_max = 1.5 * (R / Ie) * delta_e_d_max;


% Eq. 12.8-17: MCE displacement
delta_MCE_S1 = 1.5 * (R / Ie) * delta_e_m_S1;
delta_MCE_S2 = 1.5 * (R / Ie) * delta_e_m_S2;
delta_MCE_V1 = 1.5 * (R / Ie) * delta_e_m_V1;
delta_MCE_V2 = 1.5 * (R / Ie) * delta_e_m_V2;
delta_MCE_mean = 1.5 * (R / Ie) * delta_e_m_mean;
delta_MCE_max = 1.5 * (R / Ie) * delta_e_m_max;

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
