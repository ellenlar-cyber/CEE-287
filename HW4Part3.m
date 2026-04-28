%% Part 3 Compare & Compute Rc for eacxh record

%% HW4 Problem 3 - Compare Numerical Rc vs Simplified Code Equations
% T=1.0s, z=5%, alpha=-0.15, 5 Loma Prieta ground motion records
%
% Three simplified equations compared:
%   1. Miranda & Akkar (2005)       - Lecture 7 slide 20
%   2. FEMA 440 / ASCE 41           - Lecture 7 slide 21
%   3. FEMA 440A improved equation  - Lecture 7 slides 23-24

clear; clc;

%% System Parameters
T     = 1.0;     % period (s)
z     = 0.05;    % damping ratio
alpha = -0.15;   % post-yield stiffness ratio (negative = degrading)

%  Rc_numerical(i) = Ce_i / Cy_min_i  for each record i
Rc_numerical = [5, 6, 7, 8, 9];  % [Parking, SLAC-1, SLAC-2, VA-1, VA-2]
Rc_avg_numerical = mean(Rc_numerical);
record_names = {'Parking','SLAC-1','SLAC-2','VA-1','VA-2'};

%%  EQUATION 1: FEMA 440 / ASCE 41 (Rmax)
t_FEMA440  = 1 + 0.15 * log(T);                   % exponent parameter
delta_ratio = 1 + 1/abs(alpha);                    % delta_d / delta_y
Rc_FEMA440  = delta_ratio + abs(alpha)^(-t_FEMA440) / 4;

fprintf('\n EQUATION 1: FEMA 440 / ASCE 41 \n');
fprintf('  t = 1 + 0.15*ln(%.1f) = %.4f\n', T, t_FEMA440);
fprintf('  delta_d/delta_y = 1 + 1/|alpha| = %.4f\n', delta_ratio);
fprintf('  Rmax = %.4f + (%.2f)^(-%.4f)/4 = %.4f\n', ...
    delta_ratio, abs(alpha), t_FEMA440, Rc_FEMA440);

%%  EQUATION 2: FEMA 440A Improved Equation
Te   = T;                          % effective period = T for our system
d    = 5;                          % no stiffness degradation
a_440A = 1 - exp(-d * Te);
Fr_Fc  = 0;                        % Fr=0 -> Fr/Fc = 0
b_440A = 1 - Fr_Fc^2;             % = 1
gamma  = alpha;                    % post-capping slope = post-yield slope

% Term 1: (delta_c/delta_y)^a = 1^a = 1  (since delta_c = delta_y)
term1_440A = 1^a_440A;

% Term 2: zero because Fr/Fc = 0
term2_440A = b_440A * (Te / (3*abs(gamma))) * Fr_Fc * 0 * sqrt(Te);
% (delta_u - delta_r)/delta_y = 0 since delta_r = delta_u

Rc_FEMA440A = term1_440A + term2_440A;

fprintf('\n EQUATION 2: FEMA 440A \n');
fprintf('  d = %d (no stiffness degradation)\n', d);
fprintf('  a = 1 - exp(-%.0f*%.1f) = %.4f\n', d, Te, a_440A);
fprintf('  b = 1 - (Fr/Fc)^2 = 1 - 0 = %.1f\n', b_440A);
fprintf('  Term 1: (delta_c/delta_y)^a = 1^%.4f = %.4f\n', a_440A, term1_440A);
fprintf('  Term 2: = 0  (because Fr = 0 -> Fr/Fc = 0)\n');
fprintf('  Rdi = %.4f\n', Rc_FEMA440A);
fprintf('  NOTE: Rdi=1 means any R>1 is potentially unsafe per FEMA 440A\n');
fprintf('        for a bilinear system with no hardening and no residual strength.\n');

%%  SUMMARY TABLE
fprintf('\n%s\n', repmat('=',1,72));
fprintf('COMPARISON TABLE  (T=%.1fs, z=%.0f%%, alpha=%.2f)\n', T, z*100, alpha);
fprintf('%s\n', repmat('=',1,72));
fprintf('%-12s  %8s\n', 'Method', 'Rc');
fprintf('%s\n', repmat('-',1,72));

for i = 1:5
    fprintf('%-12s  %8.3f  (numerical, %s)\n', ...
        'Numerical', Rc_numerical(i), record_names{i});
end
fprintf('%-12s  %8.3f  (average of 5 records)\n', 'Avg Numer.', Rc_avg_numerical);
fprintf('%s\n', repmat('-',1,72));
fprintf('%-12s  %8.3f  (FEMA 440 / ASCE 41)\n',  'FEMA440',    Rc_FEMA440);
fprintf('%-12s  %8.3f  (FEMA 440A)\n',            'FEMA440A',   Rc_FEMA440A);
fprintf('%s\n', repmat('=',1,72));

%%  BAR CHART COMPARISON
figure('Name','Problem 3 - Rc Comparison','Position',[100 100 800 500]);

% Bars for each numerical record + average
all_Rc     = [Rc_numerical, Rc_avg_numerical, Rc_FEMA440, Rc_FEMA440A];
all_labels = [record_names, {'Avg Num.'}, {'FEMA 440'}, {'FEMA 440A'}];
colors     = [repmat([0.2 0.5 0.8], 5, 1);   % blue  - numerical records
              0.1 0.3 0.6;                    % dark blue - average
              0.2 0.7 0.3;                    % green  - FEMA 440
              0.7 0.2 0.7];                   % purple - FEMA 440A

b = bar(all_Rc, 'FaceColor','flat');
b.CData = colors;
set(gca, 'XTickLabel', all_labels, 'XTick', 1:length(all_Rc));
ylabel('R_c  (Maximum Strength Reduction Factor)');
title(sprintf('Comparison of R_c: Numerical vs Code Equations\n T = %.1fs,  \\xi = %.0f%%,  \\alpha = %.2f', ...
    T, z*100, alpha));
grid on; box on;
yline(Rc_FEMA440,      'g--', 'LineWidth', 1.5, 'Label', 'FEMA 440');
yline(1,               'm--', 'LineWidth', 1.5, 'Label', 'FEMA 440A');