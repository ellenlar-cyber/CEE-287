%% Part 1 - Plot Multi-Period DBE and MCE Spectra (0 to 2s)

% Load CSV files
dbe = readmatrix('Multi-Period Design Spectrum-20260420.csv');
mce = readmatrix('Multi-Period MCER Spectrum-20260420.csv');

% Extract columns
T_dbe = dbe(:,1);   Sa_dbe = dbe(:,2);
T_mce = mce(:,1);   Sa_mce = mce(:,2);

% Trim to 0-2s
idx_dbe = T_dbe <= 2.0;
idx_mce = T_mce <= 2.0;

% Plot
figure;
plot(T_dbe(idx_dbe), Sa_dbe(idx_dbe), 'b-', 'LineWidth', 2); hold on;
plot(T_mce(idx_mce), Sa_mce(idx_mce), 'r-', 'LineWidth', 2);
xlabel('Period T (s)');
ylabel('Spectral Acceleration S_a (g)');
title('ASCE 7-22 Multi-Period Response Spectra — Stanford Campus');
legend('DBE (Design Spectrum)', 'MCE_R Spectrum', 'Location', 'northeast');
grid on;
xlim([0 2]);
ylim([0 2.5]);

%% Part 2 - Code Displacement Estimates (Eq. 12.8-16 and 12.8-17)

% BRBF parameters from ASCE 7-22 Table 12.2-1
R  = 8;
Cd = 5;
Ie = 1.0;   % Risk Category II
T  = 0.6;   % fundamental period (s)
g  = 981;   % cm/s^2

% Interpolate Sa at T=0.6s from multi-period DBE spectrum
Sa_DBE_06 = interp1(T_dbe, Sa_dbe, T);   % in g
fprintf('Sa_DBE at T=0.6s = %.4f g\n', Sa_DBE_06);

% Interpolate Sa at T=0.6s from multi-period DBE spectrum
Sa_MCE_06 = interp1(T_mce, Sa_mce, T);   % in g
fprintf('Sa_MCE at T=0.6s = %.4f g\n', Sa_MCE_06);

% Elastic spectral displacement at DBE level (from ELF)
% delta_e = (Sa/(R/Ie)) * g*T^2/(4*pi^2)
delta_e = (Sa_DBE_06 / (R/Ie)) * g * T^2 / (4*pi^2);
fprintf('delta_e (elastic ELF displacement) = %.4f cm\n', delta_e);

% Eq. 12.8-16: DBE displacement
delta_DE = (Cd / Ie) * delta_e;
fprintf('delta_DE (Eq. 12.8-16, DBE) = %.4f cm\n', delta_DE);

% Eq. 12.8-17: MCE displacement
delta_MCE = 1.5 * (R / Ie) * delta_e;
fprintf('delta_MCE (Eq. 12.8-17, MCE) = %.4f cm\n', delta_MCE);