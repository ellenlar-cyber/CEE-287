%% Part 4 - Ground Motion Scaling to MCE Level
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

%% Load pre-computed RotD100 max-direction spectra
slac_rot       = load('SLAC_RotD100.mat');
va_rot         = load('VA_RotD100.mat');
Te_SLAC        = slac_rot.Te;
Sa_maxdir_SLAC = slac_rot.SaRotD100;
Te_VA          = va_rot.Te;
Sa_maxdir_VA   = va_rot.SaRotD100;

%% Load MCE spectrum from Part 1
mce    = readmatrix('Multi-Period MCER Spectrum-20260420.csv');
T_mce  = mce(:,1);
Sa_mce = mce(:,2);

%% Find MCE scale factors for each pair
% Period range [0.2T, 2.0T] for T = 0.6s
T_struct = 0.6;
T_lo     = 0.2 * T_struct;   % 0.12 s
T_hi     = 2.0 * T_struct;   % 1.20 s

T_common = (T_lo : 0.001 : T_hi)';

% Interpolate everything onto that grid
Sa_SLAC_common = interp1(Te_SLAC, Sa_maxdir_SLAC, T_common, 'linear', 'extrap');
Sa_VA_common   = interp1(Te_VA,   Sa_maxdir_VA,   T_common, 'linear', 'extrap');
Sa_mce_common  = interp1(T_mce,   Sa_mce,         T_common, 'linear', 'extrap');

% find average of Sa arrays
SaAvg_SLAC = mean(Sa_SLAC_common);
SaAvg_VA   = mean(Sa_VA_common);
SaAvg_mce = mean(Sa_mce_common);

% Scale SLAC and VA vs USGS
SF1 = SaAvg_mce / SaAvg_SLAC;
SF2 = SaAvg_mce / SaAvg_VA;

% Find mean for each period
Sa_mean = 0.5 * (SF1 * Sa_SLAC_common + SF2 * Sa_VA_common);

ratios = (0.9 * Sa_mce_common)./ Sa_mean;
SF3 = max(ratios);  % worst-case period sets SF3

% Final scale factors
SF_SLAC_MCE = SF1 * SF3;
SF_VA_MCE   = SF2 * SF3;

fprintf('=== MCE Scale Factors ===\n');
fprintf('SLAC pair: %.4f\n', SF_SLAC_MCE);
fprintf('VA pair:   %.4f\n', SF_VA_MCE);

%% Scale ground motion records
ag_SLAC1_sc = ag_SLAC1 * SF_SLAC_MCE;
ag_SLAC2_sc = ag_SLAC2 * SF_SLAC_MCE;
ag_VA1_sc   = ag_VA1   * SF_VA_MCE;
ag_VA2_sc   = ag_VA2   * SF_VA_MCE;

%% Compute individual response spectra of scaled records (5% damping)
T_vec = (0.05:0.05:2.0)';
z  = 0.05;

Sa_S1 = zeros(length(T_vec),1);
Sa_S2 = zeros(length(T_vec),1);
Sa_V1 = zeros(length(T_vec),1);
Sa_V2 = zeros(length(T_vec),1);

fprintf('Computing response spectra for scaled records...\n');

for i = 1:length(T_vec)
    T  = T_vec(i);
    wn = 2*pi/T;

    [~,~,~,Sd,~,~,~,~] = SDOF_Response(T, z, ag_SLAC1_sc, dt, 0, 0);
    Sa_S1(i) = wn^2 * Sd / g;

    [~,~,~,Sd,~,~,~,~] = SDOF_Response(T, z, ag_SLAC2_sc, dt, 0, 0);
    Sa_S2(i) = wn^2 * Sd / g;

    [~,~,~,Sd,~,~,~,~] = SDOF_Response(T, z, ag_VA1_sc, dt, 0, 0);
    Sa_V1(i) = wn^2 * Sd / g;

    [~,~,~,Sd,~,~,~,~] = SDOF_Response(T, z, ag_VA2_sc, dt, 0, 0);
    Sa_V2(i) = wn^2 * Sd / g;
end

Sa_mean_MCE = (Sa_S1 + Sa_S2 + Sa_V1 + Sa_V2) / 4;









%% FROM HERE IDK
% %% Print Sa values at T = 0.6s
% Sa_S1_06  = interp1(T_vec, Sa_S1,        T_struct);
% Sa_S2_06  = interp1(T_vec, Sa_S2,        T_struct);
% Sa_V1_06  = interp1(T_vec, Sa_V1,        T_struct);
% Sa_V2_06  = interp1(T_vec, Sa_V2,        T_struct);
% Sa_mean06 = interp1(T_vec, Sa_mean_MCE,  T_struct);
% 
% fprintf('\n=== Sa at T=0.6s (MCE scaled, 5%% damping) ===\n');
% fprintf('SLAC-1: %.4f g\n', Sa_S1_06);
% fprintf('SLAC-2: %.4f g\n', Sa_S2_06);
% fprintf('VA-1:   %.4f g\n', Sa_V1_06);
% fprintf('VA-2:   %.4f g\n', Sa_V2_06);
% fprintf('Mean:   %.4f g\n', Sa_mean06);

%% Plot MCE scaled spectra
idx_plot = T_vec <= 2.0;
idx_mce  = T_mce <= 2.0;

figure;
plot(T_vec(idx_plot), Sa_S1(idx_plot),       'b',   'LineWidth', 1.2); hold on;
plot(T_vec(idx_plot), Sa_S2(idx_plot),       'b--', 'LineWidth', 1.2);
plot(T_vec(idx_plot), Sa_V1(idx_plot),       'r',   'LineWidth', 1.2);
plot(T_vec(idx_plot), Sa_V2(idx_plot),       'r--', 'LineWidth', 1.2);
plot(T_mce(idx_mce),  Sa_mce(idx_mce),       'k',   'LineWidth', 2.5);
plot(T_vec(idx_plot), Sa_mean_MCE(idx_plot), 'g',   'LineWidth', 2.5);
xline(T_lo, '--k', '0.2T = 0.12s', 'LabelVerticalAlignment', 'bottom');
xline(T_hi, '--k', '2.0T = 1.20s', 'LabelVerticalAlignment', 'bottom');
xline(T_struct, '--m', 'T = 0.6s', 'LabelVerticalAlignment', 'bottom');
xlabel('Period T (s)');
ylabel('S_a (g)');
title('Scaled Ground Motion Spectra — MCE Level');
legend('SLAC-1','SLAC-2','VA-1','VA-2','MCE Target','Mean of 4 Records', ...
       'Location','northeast');
grid on;
xlim([0 2]);
ylim([0 inf]);


