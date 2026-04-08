%% Verification Script for Nonlinear SDOF
% Uses first 10 seconds of El Centro record (dt = 0.002s)
% Compares results for two systems tested in Chopra Chapter 7
% and verifies hysteretic behavior looks correct
clear; clc; close all;

%% Load El Centro Record
% Data has two columns: time (s) and acceleration (g)
data = readmatrix('First_10_secs_of_El_Centro_at_002s.xls');
t_record = data(:, 1);
ag_g     = data(:, 2);       % acceleration in g

% Convert to cm/s^2
g   = 981;                   % cm/s^2
ag  = ag_g * g;              % now in cm/s^2
dt  = t_record(2) - t_record(1);   % 0.002 s

fprintf('Record loaded: %d points, dt = %.4f s, duration = %.2f s\n', ...
    length(ag), dt, t_record(end));
fprintf('PGA = %.4f g (%.2f cm/s^2)\n', max(abs(ag_g)), max(abs(ag)));

%% Test Case 1: Tn = 0.5s, z = 5%, Cy = 0.25
% This matches examples in Chopra Chapter 7
Tn1  = 0.5;
z1   = 0.05;
Cy1  = 0.25;

[u1, ud1, udd1, Fs1, Sd1, mu1] = SDOF_Response_NL(...
    Tn1, z1, ag, dt, 0, 0, Cy1, 'Cy', 'average');

% Rebuild time vector (may have been interpolated inside function)
t1 = (0:length(u1)-1)' * (Tn1/40);

fprintf('\n--- Test Case 1: Tn=%.1fs, z=%.0f%%, Cy=%.2f ---\n', Tn1, z1*100, Cy1);
fprintf('Peak displacement  = %.4f cm\n', Sd1);
fprintf('Ductility demand   = %.4f\n',    mu1);
fprintf('Peak restoring Fs  = %.4f (normalized by W: %.4f)\n', ...
    max(abs(Fs1)), max(abs(Fs1))/g);

%% Test Case 2: Tn = 1.0s, z = 5%, Cy = 0.15
Tn2  = 1.0;
z2   = 0.05;
Cy2  = 0.15;

[u2, ud2, udd2, Fs2, Sd2, mu2] = SDOF_Response_NL(...
    Tn2, z2, ag, dt, 0, 0, Cy2, 'Cy', 'average');

t2 = (0:length(u2)-1)' * (Tn2/40);

fprintf('\n--- Test Case 2: Tn=%.1fs, z=%.0f%%, Cy=%.2f ---\n', Tn2, z2*100, Cy2);
fprintf('Peak displacement  = %.4f cm\n', Sd2);
fprintf('Ductility demand   = %.4f\n',    mu2);
fprintf('Peak restoring Fs  = %.4f (normalized by W: %.4f)\n', ...
    max(abs(Fs2)), max(abs(Fs2))/g);

%% Elastic comparison: same systems with very high Cy (elastic behavior)
% If Cy >> PGA the system should remain elastic
% In that case ductility demand should be ~1.0
Cy_elastic = 10.0;   % much larger than PGA, system stays elastic

[u1e, ~, ~, ~, Sd1e, mu1e] = SDOF_Response_NL(...
    Tn1, z1, ag, dt, 0, 0, Cy_elastic, 'Cy', 'average');
[u2e, ~, ~, ~, Sd2e, mu2e] = SDOF_Response_NL(...
    Tn2, z2, ag, dt, 0, 0, Cy_elastic, 'Cy', 'average');

fprintf('\n--- Elastic check (Cy=%.1f, should give mu~1) ---\n', Cy_elastic);
fprintf('Tn=%.1fs: mu = %.4f  (expected ~1.0)\n', Tn1, mu1e);
fprintf('Tn=%.1fs: mu = %.4f  (expected ~1.0)\n', Tn2, mu2e);

%% Also run with uy input type to verify both input paths give same result
wn1 = 2*pi/Tn1;
k1  = wn1^2;            % unit mass
Fy1 = Cy1 * 1 * g;
uy1 = Fy1 / k1;

[u1_uy, ~, ~, ~, Sd1_uy, mu1_uy] = SDOF_Response_NL(...
    Tn1, z1, ag, dt, 0, 0, uy1, 'uy', 'average');

fprintf('\n--- Input type check: Cy vs uy should give identical results ---\n');
fprintf('Sd via Cy = %.6f cm\n', Sd1);
fprintf('Sd via uy = %.6f cm\n', Sd1_uy);
fprintf('Match: %s\n', string(abs(Sd1 - Sd1_uy) < 1e-8));

%% Plotting

% Figure 1: Test Case 1 - Time Histories
figure('Name', 'Test Case 1: Tn=0.5s', 'Position', [100 100 900 700])

subplot(3,1,1)
plot(t_record, ag_g, 'k', 'LineWidth', 0.8)
xlabel('Time (s)'); ylabel('a_g (g)')
title('Ground Motion - First 10s El Centro NS')
grid on

subplot(3,1,2)
plot(t1, u1, 'b', 'LineWidth', 1)
xlabel('Time (s)'); ylabel('u (cm)')
title(sprintf('Relative Displacement  |  Tn=%.1fs, \\zeta=%.0f%%, Cy=%.2f', Tn1, z1*100, Cy1))
grid on

subplot(3,1,3)
plot(t1, Fs1/g, 'r', 'LineWidth', 1)
yline(Cy1,  '--k', sprintf('  +Cy = %.2f', Cy1))
yline(-Cy1, '--k', sprintf('  -Cy = %.2f', Cy1))
xlabel('Time (s)'); ylabel('Fs/W')
title('Normalized Restoring Force Time History')
grid on

% Figure 2: Test Case 1 - Hysteretic Behavior
figure('Name', 'Hysteretic Behavior', 'Position', [100 100 900 400])

subplot(1,2,1)
plot(u1, Fs1/g, 'b', 'LineWidth', 1)
xlabel('Displacement u (cm)'); ylabel('Fs/W')
title(sprintf('Tn=%.1fs, Cy=%.2f,  \\mu=%.2f', Tn1, Cy1, mu1))
grid on

subplot(1,2,2)
plot(u2, Fs2/g, 'r', 'LineWidth', 1)
xlabel('Displacement u (cm)'); ylabel('Fs/W')
title(sprintf('Tn=%.1fs, Cy=%.2f,  \\mu=%.2f', Tn2, Cy2, mu2))
grid on

% Figure 3: Test Case 2 - Time Histories
figure('Name', 'Test Case 2: Tn=1.0s', 'Position', [100 100 900 500])

subplot(2,1,1)
plot(t2, u2, 'r', 'LineWidth', 1)
xlabel('Time (s)'); ylabel('u (cm)')
title(sprintf('Relative Displacement  |  Tn=%.1fs, \\zeta=%.0f%%, Cy=%.2f', Tn2, z2*100, Cy2))
grid on

subplot(2,1,2)
plot(t2, Fs2/g, 'r', 'LineWidth', 1)
yline(Cy2,  '--k', sprintf('  +Cy = %.2f', Cy2))
yline(-Cy2, '--k', sprintf('  -Cy = %.2f', Cy2))
xlabel('Time (s)'); ylabel('Fs/W')
title('Normalized Restoring Force Time History')
grid on

% Figure 4: Elastic vs Inelastic displacement comparison
figure('Name', 'Elastic vs Inelastic', 'Position', [100 100 900 400])

t1e = (0:length(u1e)-1)' * (Tn1/40);
t2e = (0:length(u2e)-1)' * (Tn2/40);

subplot(1,2,1)
plot(t1e, u1e, 'k--', 'LineWidth', 1); hold on
plot(t1,  u1,  'b',   'LineWidth', 1)
legend('Elastic (Cy=10)', sprintf('Inelastic (Cy=%.2f)', Cy1))
xlabel('Time (s)'); ylabel('u (cm)')
title(sprintf('Tn=%.1fs', Tn1))
grid on

subplot(1,2,2)
plot(t2e, u2e, 'k--', 'LineWidth', 1); hold on
plot(t2,  u2,  'r',   'LineWidth', 1)
legend('Elastic (Cy=10)', sprintf('Inelastic (Cy=%.2f)', Cy2))
xlabel('Time (s)'); ylabel('u (cm)')
title(sprintf('Tn=%.1fs', Tn2))
grid on

fprintf('\nAll plots generated. Check hysteretic loops for correct EPP shape.\n')
fprintf('The loops should show:\n')
fprintf('  - Linear loading/unloading slopes equal to k\n')
fprintf('  - Flat plateaus at +Fy and -Fy during yielding\n')
fprintf('  - Shifted elastic range (2uy wide) after each unloading event\n')