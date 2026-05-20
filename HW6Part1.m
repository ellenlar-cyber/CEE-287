%% HW6 Part 1 — Dynamic Properties of 9-Story Shear Building
clc; clear; close all;

%% Building Parameters
N      = 9;           % number of stories
W_tot  = 16200;       % total building weight (kips)
k_story = 1700;       % uniform story stiffness (kip/in)
g      = 386.1;       % gravitational acceleration (in/s^2)
zeta   = 0.035;       % damping ratio (used in later parts)

% Uniform floor mass (kip*s^2/in) — same at every floor
m_floor = W_tot / (N * g);

%% Step 1 Build mass and stiffness vectors
% delta=1 gives uniform distribution; lambda irrelevant when delta=1
m_vec = mass_vector(N, 1, 1) * m_floor;   % (9x1) uniform mass
k_vec = stiffness_vector(N, 1, 1) * k_story; % (9x1) uniform stiffness

%% Step 2 Run eigenvalue analysis via mdof_analysis
% Returns R.T (all 9 periods), R.PHI (normalized modes 1-3 only), R.M, R.K
R = mdof_analysis(N, m_vec, k_vec, 'HW6 Part 1: 9-Story Shear Building');

%% Step 3 Extend to all 9 modes: normalize PHI and compute Gamma, EMS
% mdof_analysis only normalizes modes 1-3 to unity at roof; do all 9 here
M   = R.M;
PHI = R.PHI;   % (9x9) — columns are mode shapes from eig, sorted by frequency

% Normalize all 9 mode shapes so that phi(roof, n) = 1
PHI_norm = zeros(N, N);
for n = 1:N
    PHI_norm(:,n) = PHI(:,n) / PHI(N,n);
end

% Modal participation factors for all 9 modes
% Gamma_n = (phi_n' * M * 1) / (phi_n' * M * phi_n)
ones_vec = ones(N,1);
Gamma_all = zeros(N,1);
for n = 1:N
    L_n          = PHI_norm(:,n)' * M * ones_vec;
    M_n          = PHI_norm(:,n)' * M * PHI_norm(:,n);
    Gamma_all(n) = L_n / M_n;
end

% Effective mode shapes for all 9 modes: EMS(:,n) = Gamma_n * phi_n
EMS_all = zeros(N, N);
for n = 1:N
    EMS_all(:,n) = Gamma_all(n) * PHI_norm(:,n);
end

%% Print table of periods and period ratios
fprintf('\n Periods and Period Ratios\n');
fprintf('%-6s %-12s %-12s\n','Mode','T_n (s)','T_n/T_1');
T1 = R.T(1);
for n = 1:N
    fprintf('%-6d %-12.4f %-12.4f\n', n, R.T(n), R.T(n)/T1);
end

fprintf('\n Modal Participation Factors (all 9 modes) \n');
for n = 1:N
    fprintf('  Gamma_%d = %.4f\n', n, Gamma_all(n));
end

%% Plot all 9 effective mode shapes
floors = (1:N)';
figure('Name','HW6 Part 1 — All 9 Effective Mode Shapes');
for n = 1:N
    subplot(3,3,n);
    plot(EMS_all(:,n), floors, 'b-o', 'LineWidth', 1.5, 'MarkerSize', 5);
    hold on;
    xline(0,'k--');
    xlabel('\Gamma_n\phi_{jn}');
    ylabel('Floor');
    title(sprintf('Mode %d  (T=%.3fs)', n, R.T(n)));
    grid on;
    ylim([0 N+0.5]);
end
sgtitle('Effective Mode Shapes — 9-Story Shear Building');

%% Save outputs for use in later parts
save('HW6_Part1_results.mat', 'R', 'PHI_norm', 'Gamma_all', 'EMS_all', ...
     'M', 'm_vec', 'k_vec', 'N', 'W_tot', 'g', 'zeta', 'm_floor', 'k_story');