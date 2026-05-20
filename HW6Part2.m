%% HW6 Part 2 — Response Spectrum + Equivalent Static Analysis (ESA)
clc; clear; close all;

%% Load Part 1 results
load('HW6_Part1_results.mat');
% Loaded: R, PHI_norm, Gamma_all, EMS_all, M, m_vec, k_vec,
%         N, W_tot, g, zeta, m_floor, k_story
floors  = (1:N)';
stories = (1:N)';

% Building height
H_total = 118 * 12;        % total height (in)
h_story = H_total / N;     % uniform story height (in) = 157.3 in

%% Step 1 Load ground motion
% readmatrix skips blank and header rows automatically
data = readmatrix('gm_hw_6.xlsx');
T_modes  = R.T;
% call the 5th modal period in order to get necessary time step
% will be Tn/20 of the 5th modal period
dt = T_modes(5)/20;
% acceleration time history (cm/s^2)
ag   = data(:,2);
% time
t_old = data(:,1);
% new time vector
t_new = 0:dt:t_old(end);  % new time vector based on the time step
ag = interp1(t_old, ag, t_new, 'linear');  % interpolate ground motion

%% Step 2  Compute response spectrum (Lecture 1 Iwan's method in SDOF_Response)
% Periods: 0.02 to 3.0s in steps of 0.02s = 150 periods
T_spec  = (0.02:0.02:3.0)';
nT      = length(T_spec);

Spa_spec = zeros(nT,1);   % pseudo-acceleration spectrum (cm/s^2)
Sd_spec  = zeros(nT,1);   % displacement spectrum (cm)

for i = 1:nT
    [~,~,~, Sd_spec(i),~,~,~, Spa_spec(i)] = SDOF_Response(T_spec(i), zeta, ag, dt, 0, 0);
end

%% Step 3 Interpolate Spa and Sd at the 9 modal periods
T_modes  = R.T;   % (9x1) modal periods from Part 1

Spa_modes = interp1(T_spec, Spa_spec, T_modes, 'linear');  % cm/s^2
Sd_modes  = interp1(T_spec, Sd_spec,  T_modes, 'linear');  % cm

% Convert to inch units for building calculations
Spa_modes_in = Spa_modes / 2.54;   % in/s^2
Sd_modes_in  = Sd_modes  / 2.54;   % inches

% Print table
fprintf('\n Spectral Ordinates at 9 Modal Periods\n');
fprintf('%-6s %-10s %-14s %-14s %-14s %-14s\n', ...
    'Mode','T_n (s)','Spa (cm/s^2)','Sd (cm)','Spa (in/s^2)','Sd (in)');
for n = 1:N
    fprintf('%-6d %-10.4f %-14.4f %-14.4f %-14.4f %-14.4f\n', ...
        n, T_modes(n), Spa_modes(n), Sd_modes(n), Spa_modes_in(n), Sd_modes_in(n));
end

%% Step 4 ESA using mode 1 only (R=1, unreduced)
phi1    = PHI_norm(:,1);        % first mode shape (normalized to 1 at roof)
Gamma1  = Gamma_all(1);         % modal participation factor, mode 1

Spa1    = Spa_modes_in(1);      % Spa at T1 (in/s^2)
Sd1     = Sd_modes_in(1);       % Sd  at T1 (inches)

% store beta 2 for mode 1
Beta2_1 = R.beta2;
Beta1_1 = R.beta1;

% Effective modal mass M1* = L1^2 / Mmodal1  (kip*s^2/in)
L1      = phi1' * M * ones(N,1);
Mmodal1 = phi1' * M * phi1;
M1star  = L1^2 / Mmodal1;

% Base shear
% V = CsW
V_ESA   = M1star * Spa1;                  % kips
VW_ESA  = V_ESA / W_tot;                  % fraction of total weight

% Lateral forces at each floor
% Fjn = gamma n * phi jn * mj *An
F_ESA   = m_vec .* (Gamma1 * phi1) * Spa1;  % kips (Nx1)

% Story shears (sum F from floor x to roof)
% sum for story forces at that level
Vs_ESA  = zeros(N,1);
for x = 1:N
    Vs_ESA(x) = sum(F_ESA(x:N));
end

% Floor displacements
u_ESA   = Gamma1 * phi1 * Sd1;            % inches (Nx1)

% Interstory drift ratios
IDR_ESA = zeros(N,1);
IDR_ESA(1) = u_ESA(1) / h_story;          % story 1: relative to ground
for x = 2:N
    IDR_ESA(x) = (u_ESA(x) - u_ESA(x-1)) / h_story;
end

% Print ESA results
fprintf('\nESA Results\n');
fprintf('Base shear V       = %.2f kips\n', V_ESA);
fprintf('V / W_total        = %.4f (%.2f%%)\n', VW_ESA, VW_ESA*100);
fprintf('\n%-6s %-12s %-14s %-14s %-14s\n', ...
    'Floor','F_i (kips)','V_story (kips)','u_i (in)','IDR_i (%)');
for i = 1:N
    fprintf('%-6d %-12.4f %-14.4f %-14.4f %-14.4f\n', ...
        i, F_ESA(i), Vs_ESA(i), u_ESA(i), IDR_ESA(i)*100);
end

%% Step 5 Plots

% Figure 1: Response Spectrum
figure('Name','HW6 Part 2: Response Spectrum');
plot(T_spec, Spa_spec/981, 'b-', 'LineWidth', 1.5);
hold on;
plot(T_modes, Spa_modes/981, 'ro', 'MarkerSize', 8, 'MarkerFaceColor','r');
for n = 1:N
    text(T_modes(n), Spa_modes(n)/981 + 0.02, sprintf('T_%d',n), ...
        'FontSize', 7, 'HorizontalAlignment','center');
end
xlabel('Period T (s)');
ylabel('S_{pa} (g)');
title('Linear Elastic Response Spectrum : \zeta = 3.5%');
legend('Spectrum','Modal periods','Location','NorthEast');
grid on;
xlim([0 3.0]);

% Figure 2: ESA Response Profiles
figure('Name','HW6 Part 2: ESA Profiles');

% Lateral forces
subplot(1,4,1);
barh(floors, F_ESA, 0.5, 'FaceColor',[0.2 0.5 0.8]);
xlabel('F_i (kips)');
ylabel('Floor');
title('Lateral Forces');
grid on; ylim([0 N+1]);

% Story shears
subplot(1,4,2);
plot(Vs_ESA, stories, 'b-o', 'LineWidth', 1.5, 'MarkerFaceColor','b');
xlabel('V_x (kips)');
ylabel('Story');
title('Story Shears');
grid on; ylim([0 N+1]);

% Floor displacements
subplot(1,4,3);
plot([0; u_ESA], [0; floors], 'b-o', 'LineWidth', 1.5, 'MarkerFaceColor','b');
xlabel('u_i (in)');
ylabel('Floor');
title('Floor Displacements');
grid on; ylim([0 N+1]);

% Interstory drift ratios
subplot(1,4,4);
plot(IDR_ESA*100, stories, 'b-o', 'LineWidth', 1.5, 'MarkerFaceColor','b');
xlabel('IDR (%)');
ylabel('Story');
title('Interstory Drift Ratios');
grid on; ylim([0 N+1]);

sgtitle(sprintf('ESA Results: Mode 1 Only, R=1  (V = %.1f kips = %.2f%% W)', ...
    V_ESA, VW_ESA*100));

%% Save outputs for Parts 3 and 4
save('HW6_Part2_results.mat', ...
     'T_spec','Spa_spec','Sd_spec', ...
     'Spa_modes','Sd_modes','Spa_modes_in','Sd_modes_in', ...
     'V_ESA','VW_ESA','F_ESA','Vs_ESA','u_ESA','IDR_ESA', ...
     'ag','dt','Beta1_1','Beta2_1');