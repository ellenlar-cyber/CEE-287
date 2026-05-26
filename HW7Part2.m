%% HW7 Part B: Peak Floor Accelerations at All Levels (4 Methods)
clc; clear; close all;

%% Load Part 1 results
load('HW7_Part1_results.mat');
% Loaded: T_approx, PHI_approx, Gamma_approx,
%         a_floor_abs, a_floor_rel, a_recorded,
%         ag, dt, t, N, H_total, h_story, x_floors, zeta, T1, Nmodes

floors = (1:N)';

%% Method 1 — RHA peak, relative SDOF accelerations
% peak is max absolute value over entire response history (reused from HW6 Part 4)
PFA_rel = max(abs(a_floor_rel), [], 1)';    % (9x1) cm/s^2

%% Method 2 — RHA peak, absolute SDOF accelerations
PFA_abs = max(abs(a_floor_abs), [], 1)';    % (9x1) cm/s^2

%% Method 3 — RSA with absolute acceleration response spectrum of base motion
% Compute 5%-damped absolute acceleration spectrum of recorded base motion
% over a fine period grid, then interpolate at the three modal periods
T_spec  = (0.02:0.02:4.0)';
nT      = length(T_spec);
Sa_abs_spec = zeros(nT, 1);

for i = 1:nT
    % udd_abs (3rd output) is the absolute acceleration — same call as Part 1
    [~,~, udd_abs_i, ~,~,~,~,~] = SDOF_Response(T_spec(i), zeta, ag, dt, 0, 0);
    Sa_abs_spec(i) = max(abs(udd_abs_i));   % peak absolute acc (cm/s^2)
end

% Interpolate at the three approximate modal periods
Sa_abs_modes = interp1(T_spec, Sa_abs_spec, T_approx, 'linear');  % (3x1) cm/s^2

% SRSS combination over 3 modes at each floor
% PFA_RSA(j) = sqrt( sum_n [ Gamma_n * phi_n(j) * Sa_abs(Tn) ]^2 )
PFA_RSA = zeros(N, 1);
for j = 1:N
    modal_contributions = zeros(Nmodes, 1);
    for n = 1:Nmodes
        modal_contributions(n) = Gamma_approx(n) * PHI_approx(j,n) * Sa_abs_modes(n);
    end
    PFA_RSA(j) = sqrt(sum(modal_contributions.^2));
end

%% Method 4 — ASCE 7-22 code distribution (Chapter 13)
% PFA(z) = PGA * (1 + 2*z/H)  — linear from 1xPGA at base to 3xPGA at roof
PGA = max(abs(ag));     % peak ground acceleration from recorded base (cm/s^2)

PFA_code = PGA * (1 + 2 * x_floors / H_total);   % (9x1) cm/s^2

%% Print results table
fprintf('HW7 Part 2 — Peak Floor Accelerations (cm/s^2)\n');
fprintf('%-6s %-14s %-14s %-14s %-14s\n', ...
    'Floor', 'RHA-Rel', 'RHA-Abs', 'RSA', 'Code');
fprintf('%s\n', repmat('-',1,62));
for j = 1:N
    fprintf('%-6d %-14.2f %-14.2f %-14.2f %-14.2f\n', ...
        j, PFA_rel(j), PFA_abs(j), PFA_RSA(j), PFA_code(j));
end
fprintf('\nPGA (base) = %.2f cm/s^2\n', PGA);
fprintf('Sa_abs at modal periods:\n');
for n = 1:Nmodes
    fprintf('  Mode %d: T = %.3f s,  Sa_abs = %.2f cm/s^2\n', ...
        n, T_approx(n), Sa_abs_modes(n));
end

%% Plot — PFA profile vs floor height
figure('Name','HW7 Part 2: Peak Floor Accelerations');

plot(PFA_rel,  floors, 'b-o',  'LineWidth', 1.5, 'MarkerFaceColor', 'b'); hold on;
plot(PFA_abs,  floors, 'r-o',  'LineWidth', 1.5, 'MarkerFaceColor', 'r');
plot(PFA_RSA,  floors, 'g-s',  'LineWidth', 1.5, 'MarkerFaceColor', 'g');
plot(PFA_code, floors, 'k--^', 'LineWidth', 1.5, 'MarkerFaceColor', 'k');
xline(PGA, 'k:', 'LineWidth', 1.0);   % mark PGA for reference

xlabel('Peak Floor Acceleration (cm/s^2)');
ylabel('Floor');
title('Peak Floor Accelerations — All 9 Levels');
legend('RHA: Relative SDOF acc', ...
       'RHA: Absolute SDOF acc', ...
       'RSA: Absolute acc spectrum', ...
       'ASCE 7-22 Code', ...
       'PGA', ...
       'Location', 'SouthEast');
grid on;
ylim([0 N+1]);
yticks(1:N);

%% Save for Part 3
save('HW7_Part2_results.mat', ...
    'PFA_rel', 'PFA_abs', 'PFA_RSA', 'PFA_code', ...
    'PGA', 'Sa_abs_spec', 'Sa_abs_modes', 'T_spec');