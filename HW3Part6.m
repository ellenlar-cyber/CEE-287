%% Part 6 Displacement demands by inelastic displacement ratios
clear; clc; close all;

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

% Scaled factors 
SF_SLAC_DBE = 2.8634;
SF_VA_DBE =  2.7373;
SF_SLAC_MCE = 4.2972;
SF_VA_MCE  =  4.1080;

%%  DBE Scale ground motion records
ag_SLAC1_sc_dbe = ag_SLAC1 * SF_SLAC_DBE;
ag_SLAC2_sc_dbe = ag_SLAC2 * SF_SLAC_DBE;
ag_VA1_sc_dbe   = ag_VA1   * SF_VA_DBE;
ag_VA2_sc_dbe   = ag_VA2   * SF_VA_DBE;

%% MCE Scale ground motion records
ag_SLAC1_sc_mce = ag_SLAC1 * SF_SLAC_MCE;
ag_SLAC2_sc_mce = ag_SLAC2 * SF_SLAC_MCE;
ag_VA1_sc_mce   = ag_VA1   * SF_VA_MCE;
ag_VA2_sc_mce   = ag_VA2   * SF_VA_MCE;

%% Parameters
T = 0.6;
z  = 0.02;
R = 8;
Cd = 5;
Ie = 1.0;
z2 = 0.02;   % 2% damping for RHA
z5 = 0.05;   % 5% damping for Sa 
wn = 2*pi/T;

a = 60; % Site class D

%% Compute Sa values for 5%
[~,~,~,~,~,Sa5,~,~] = SDOF_Response(T, z5, ag_SLAC1_sc_dbe, dt, 0, 0);
Sa_S1_dbe = Sa5/g;

[~,~,~,~,~,Sa5,~,~] = SDOF_Response(T, z5, ag_SLAC2_sc_dbe, dt, 0, 0);
Sa_S2_dbe = Sa5/g;


[~,~,~,~,~,Sa5,~,~] = SDOF_Response(T, z5, ag_VA1_sc_dbe, dt, 0, 0);
Sa_V1_dbe = Sa5/g;

[~,~,~,~,~,Sa5,~,~] = SDOF_Response(T, z5, ag_VA2_sc_dbe, dt, 0, 0);
Sa_V2_dbe = Sa5/g;


% mu_strength per record
% mu is Sa(of record)/(Vy/W)
% Sa is taken from 5% damping case, calculated above
% Vy/W = Sa design /R
Sa_design_dbe = 1.4800;
Vy_dbe = Sa_design_dbe / R;

% Cm = 0 because "other" structure
Cm = 1.0;

mu_S1_dbe = Sa_S1_dbe / Vy_dbe * Cm;
mu_S2_dbe = Sa_S2_dbe / Vy_dbe * Cm;
mu_V1_dbe = Sa_V1_dbe / Vy_dbe * Cm;
mu_V2_dbe = Sa_V2_dbe / Vy_dbe * Cm;

% C1 per record
C1_S1_dbe = 1 + (mu_S1_dbe - 1) / (a * T^2);
C1_S2_dbe = 1 + (mu_S2_dbe - 1) / (a * T^2);
C1_V1_dbe = 1 + (mu_V1_dbe - 1) / (a * T^2);
C1_V2_dbe = 1 + (mu_V2_dbe - 1) / (a * T^2);

%% Sd at 2% damping, multiply by C1
[~,~,~,Sd2,~,~,~,~] = SDOF_Response(T, z2, ag_SLAC1_sc_dbe, dt, 0, 0);
Sd_S1_dbe   = Sd2;
Disp_S1_dbe = C1_S1_dbe * Sd2;

[~,~,~,Sd2,~,~,~,~] = SDOF_Response(T, z2, ag_SLAC2_sc_dbe, dt, 0, 0);
Sd_S2_dbe   = Sd2;
Disp_S2_dbe = C1_S2_dbe * Sd2;

[~,~,~,Sd2,~,~,~,~] = SDOF_Response(T, z2, ag_VA1_sc_dbe, dt, 0, 0);
Sd_V1_dbe   = Sd2;
Disp_V1_dbe = C1_V1_dbe * Sd2;

[~,~,~,Sd2,~,~,~,~] = SDOF_Response(T, z2, ag_VA2_sc_dbe, dt, 0, 0);
Sd_V2_dbe   = Sd2;
Disp_V2_dbe = C1_V2_dbe * Sd2;

Disp_mean_dbe = (Disp_S1_dbe + Disp_S2_dbe + Disp_V1_dbe + Disp_V2_dbe) / 4;
Disp_max_dbe  = max([Disp_S1_dbe, Disp_S2_dbe, Disp_V1_dbe, Disp_V2_dbe]);

% %% MCE
% [~,~,~,Sd5,~,~,~,~] = SDOF_Response(T, z5, ag_SLAC1_sc_mce, dt, 0, 0);
% Sa_S1_mce = wn^2 * Sd5 / g;
% 
% [~,~,~,Sd5,~,~,~,~] = SDOF_Response(T, z5, ag_SLAC2_sc_mce, dt, 0, 0);
% Sa_S2_mce = wn^2 * Sd5 / g;
% 
% [~,~,~,Sd5,~,~,~,~] = SDOF_Response(T, z5, ag_VA1_sc_mce, dt, 0, 0);
% Sa_V1_mce = wn^2 * Sd5 / g;
% 
% [~,~,~,Sd5,~,~,~,~] = SDOF_Response(T, z5, ag_VA2_sc_mce, dt, 0, 0);
% Sa_V2_mce = wn^2 * Sd5 / g;
% 
% % mu_strength per record
% mu_S1_mce = Sa_S1_mce / lat_strength;
% mu_S2_mce = Sa_S2_mce / lat_strength;
% mu_V1_mce = Sa_V1_mce / lat_strength;
% mu_V2_mce = Sa_V2_mce / lat_strength;
% 
% % C1 per record
% C1_S1_mce = 1 + (mu_S1_mce - 1) / (a * T^2);
% C1_S2_mce = 1 + (mu_S2_mce - 1) / (a * T^2);
% C1_V1_mce = 1 + (mu_V1_mce - 1) / (a * T^2);
% C1_V2_mce = 1 + (mu_V2_mce - 1) / (a * T^2);
% 
% %% Sd at 2% damping, multiply by C1
% 
% [~,~,~,Sd2,~,~,~,~] = SDOF_Response(T, z2, ag_SLAC1_sc_mce, dt, 0, 0);
% Sd_S1_mce   = Sd2;
% Disp_S1_mce = C1_S1_mce * Sd2;
% 
% [~,~,~,Sd2,~,~,~,~] = SDOF_Response(T, z2, ag_SLAC2_sc_mce, dt, 0, 0);
% Sd_S2_mce   = Sd2;
% Disp_S2_mce = C1_S2_mce * Sd2;
% 
% [~,~,~,Sd2,~,~,~,~] = SDOF_Response(T, z2, ag_VA1_sc_mce, dt, 0, 0);
% Sd_V1_mce   = Sd2;
% Disp_V1_mce = C1_V1_mce * Sd2;
% 
% [~,~,~,Sd2,~,~,~,~] = SDOF_Response(T, z2, ag_VA2_sc_mce, dt, 0, 0);
% Sd_V2_mce   = Sd2;
% Disp_V2_mce = C1_V2_mce * Sd2;
% 
% Disp_mean_mce = (Disp_S1_mce + Disp_S2_mce + Disp_V1_mce + Disp_V2_mce) / 4;
% Disp_max_mce  = max([Disp_S1_mce, Disp_S2_mce, Disp_V1_mce, Disp_V2_mce]);

%% Print summary
fprintf('Part 6 Results: C1 Method (ASCE 41) \n\n');
fprintf('%-12s %8s %8s %8s %8s\n', '', 'SLAC-1', 'SLAC-2', 'VA-1', 'VA-2');
fprintf('%-12s %8.4f %8.4f %8.4f %8.4f\n', 'mu_str DBE', mu_S1_dbe, mu_S2_dbe, mu_V1_dbe, mu_V2_dbe);
fprintf('%-12s %8.4f %8.4f %8.4f %8.4f\n', 'C1 DBE',     C1_S1_dbe, C1_S2_dbe, C1_V1_dbe, C1_V2_dbe);
fprintf('%-12s %8.4f %8.4f %8.4f %8.4f\n', 'Disp DBE',   Disp_S1_dbe, Disp_S2_dbe, Disp_V1_dbe, Disp_V2_dbe);
fprintf('Mean DBE = %.4f cm,  Max DBE = %.4f cm\n\n', Disp_mean_dbe, Disp_max_dbe);

% fprintf('%-12s %8.4f %8.4f %8.4f %8.4f\n', 'mu_str MCE', mu_S1_mce, mu_S2_mce, mu_V1_mce, mu_V2_mce);
% fprintf('%-12s %8.4f %8.4f %8.4f %8.4f\n', 'C1 MCE',     C1_S1_mce, C1_S2_mce, C1_V1_mce, C1_V2_mce);
% fprintf('%-12s %8.4f %8.4f %8.4f %8.4f\n', 'Disp MCE',   Disp_S1_mce, Disp_S2_mce, Disp_V1_mce, Disp_V2_mce);
% fprintf('Mean MCE = %.4f cm,  Max MCE = %.4f cm\n\n', Disp_mean_mce, Disp_max_mce);

% fprintf('Peak Elastic Displacement Demand (2%% damping)\n\n');
% fprintf('%-12s %8s %8s %8s %8s %8s %8s\n','','SLAC-1','SLAC-2','VA-1','VA-2','Mean','Max');
% fprintf('%-12s %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f\n', 'Sd DBE (cm)', ...
%     Sd_S1_dbe, Sd_S2_dbe, Sd_V1_dbe, Sd_V2_dbe, ...
%     mean([Sd_S1_dbe, Sd_S2_dbe, Sd_V1_dbe, Sd_V2_dbe]), ...
%     max([Sd_S1_dbe, Sd_S2_dbe, Sd_V1_dbe, Sd_V2_dbe]));
% % fprintf('%-12s %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f\n', 'Sd MCE (cm)', ...
% %     Sd_S1_mce, Sd_S2_mce, Sd_V1_mce, Sd_V2_mce, ...
% %     mean([Sd_S1_mce, Sd_S2_mce, Sd_V1_mce, Sd_V2_mce]), ...
% %     max([Sd_S1_mce, Sd_S2_mce, Sd_V1_mce, Sd_V2_mce]));