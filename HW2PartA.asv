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
z1   = 0.00;
Cy1  = 0.1710;

[u1, ud1, udd_abs1, Fs1, Sd1, mu1] = SDOF_Response_NL_1(...
    Tn1, z1, ag, dt, 0, 0, Cy1, 'Cy', 'linear');

% Rebuild time vector (may have been interpolated inside function)
t1 = (0:length(u1)-1)' * dt;

fprintf('\n--- Test Case 1: Tn=%.1fs, z=%.0f%%, Cy=%.2f ---\n', Tn1, z1*100, Cy1);
fprintf('Peak displacement  = %.4f cm\n', Sd1);
fprintf('Ductility demand   = %.4f\n',    mu1);
fprintf('Peak restoring Fs  = %.4f (normalized by W: %.4f)\n', ...
    max(abs(Fs1)), max(abs(Fs1))/g);


%% Also run with uy input type to verify both input paths give same result
wn1 = 2*pi/Tn1;
k1  = wn1^2;            % unit mass
Fy1 = Cy1 * 1 * g;
uy1 = Fy1 / k1;

[u1_uy, ~, ~, ~, Sd1_uy, mu1_uy] = SDOF_Response_NL_1(...
    Tn1, z1, ag, dt, 0, 0, uy1, 'uy', 'average');

fprintf('\n--- Input type check: Cy vs uy should give identical results ---\n');
fprintf('Sd via Cy = %.6f cm\n', Sd1);
fprintf('Sd via uy = %.6f cm\n', Sd1_uy);
fprintf('Match: %s\n', string(abs(Sd1 - Sd1_uy) < 1e-8));

%% Plotting

% Figure 1: Test Case 1 - Time Histories
figure('Name', 'Test Case 1: Tn=0.5s', 'Position', [100 100 900 700])

subplot(5,1,1)
plot(t_record, ag_g, 'k', 'LineWidth', 0.8)
xlabel('Time (s)'); ylabel('a_g (g)')
title('Ground Motion - First 10s El Centro NS')
grid on

subplot(5,1,2)
plot(t1, u1, 'b', 'LineWidth', 1)
xlabel('Time (s)'); ylabel('u (cm)')
title(sprintf('Relative Displacement  |  Tn=%.1fs, \\zeta=%.0f%%, Cy=%.2f', Tn1, z1*100, Cy1))
grid on

subplot(5,1,3)
plot(t1, ud1, 'r', 'LineWidth', 1)
xlabel('Time (s)'); ylabel('cm/s')
title('Relative Velocity Time History')
grid on

subplot(5,1,4)
plot(t1, udd_abs1, 'r', 'LineWidth', 1)
xlabel('Time (s)'); ylabel('cm/s^2')
title('Absolute Acceleration Time History')
grid on

% Figure 2: Test Case 1 - Hysteretic Behavior
figure('Name', 'Hysteretic Behavior', 'Position', [100 100 900 400])


plot(u1, Fs1/g, 'b', 'LineWidth', 1)
xlabel('Displacement u (cm)'); ylabel('Fs/W')
title(sprintf('Tn=%.1fs, Cy=%.2f,  \\mu=%.2f', Tn1, Cy1, mu1))
grid on


fprintf('Peak displacement:  MATLAB = %.4f cm (%.4f in)\n', Sd1, Sd1/2.54)
fprintf('Chopra Fig 7.4.2:   um = 1.71 in = 4.343 cm\n')
fprintf('Difference: %.2f%%\n', 100*abs(Sd1 - 4.343)/4.343)

%% Figure 7.4.2c - Yielding Intervals
tol = 0.001;  % tolerance for detecting yield plateau

yield_pos = abs(Fs1/g - Cy1) < tol;   % FLAG = +1
yield_neg = abs(Fs1/g + Cy1) < tol;   % FLAG = -1

figure('Name', 'Yielding Intervals', 'Position', [100 100 900 300])

subplot(2,1,1)
area(t1, yield_pos, 'FaceColor', 'k', 'EdgeColor', 'none')
ylabel('+Yield')
ylim([-0.5 1.5])
yticks([])
title('Time Intervals of Yielding')
grid on

subplot(2,1,2)
area(t1, yield_neg, 'FaceColor', 'k', 'EdgeColor', 'none')
ylabel('-Yield')
ylim([-0.5 1.5])
yticks([])
xlabel('Time (s)')
grid on

