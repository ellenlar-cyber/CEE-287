%% HW 3 Part 5 - Compute Displacement Demands

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

%% Load pre-computed RotD100 max-direction spectra
slac_rot       = load('SLAC_RotD100.mat');
va_rot         = load('VA_RotD100.mat');
Te_SLAC        = slac_rot.Te;
Sa_maxdir_SLAC = slac_rot.SaRotD100;
Te_VA          = va_rot.Te;
Sa_maxdir_VA   = va_rot.SaRotD100;

%% Load DBE spectrum from Part 1
dbe    = readmatrix('Multi-Period Design Spectrum-20260420.csv');
T_dbe  = dbe(:,1);
Sa_dbe = dbe(:,2);

%% Final scale factors DBE
SF_SLAC_DBE = 2.8634;
SF_VA_DBE   = 2.7373;

%% Final scale factors MCE
SF_SLAC_MCE = 4.2972;
SF_VA_MCE   = 4.1080;





% Eq. 12.8-16: DBE displacement
delta_DE = (Cd / Ie) * delta_e;
fprintf('delta_DE (Eq. 12.8-16, DBE) = %.4f cm\n', delta_DE);

% Eq. 12.8-17: MCE displacement
delta_MCE = 1.5 * (R / Ie) * delta_e;
fprintf('delta_MCE (Eq. 12.8-17, MCE) = %.4f cm\n', delta_MCE);