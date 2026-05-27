%% HW7 Part A: Approximate Floor Acceleration Time Histories
clc; clear; close all;

%% Building parameters (same as HW6)
N       = 9;
H_total = 118 * 12;          % total height (in)
h_story = H_total / N;       % 157.33 in per story
zeta    = 0.035;              % damping ratio (all modes)
T1      = 2.0;                % fundamental period (s) — given
Nmodes  = 3;                  % first three modes only

%% Load ground motion from hw7_acc.csv
% Columns: TIME [s], Base, 3rd floor, 7th floor, Roof  (cm/s^2)
data       = readmatrix('hw7_acc.csv', 'NumHeaderLines', 3);
t          = data(:,1);
ag         = data(:,2);         % base acceleration (cm/s^2) — SDOF input
a_3rd_rec  = data(:,3);         % recorded 3rd floor (cm/s^2)
a_7th_rec  = data(:,4);         % recorded 7th floor (cm/s^2)
a_roof_rec = data(:,5);         % recorded roof      (cm/s^2)
dt         = t(2) - t(1);       % 0.02 s

%% Approximate dynamic characteristics — uniform shear beam
% Period ratios: Ti/T1 = 1/(2i-1)
T_approx = T1 ./ (2*(1:Nmodes)' - 1);   % [2.0; 0.667; 0.400] s

% Floor heights
x_floors = (1:N)' * h_story;            % (9x1) height at each floor (in)

% Mode shapes: phi_i(x) = (-1)^(i-1) * sin( (2i-1)*pi*x / (2H) )
PHI_approx = zeros(N, Nmodes);
for i = 1:Nmodes
    PHI_approx(:,i) = (-1)^(i-1) * sin((2*i-1)*pi*x_floors / (2*H_total));
end

% Modal participation factors: Gamma_i = 4*(-1)^(i-1) / (2*i*pi - pi)
Gamma_approx = zeros(Nmodes, 1);
for i = 1:Nmodes
    Gamma_approx(i) = 4*(-1)^(i-1) / (2*i*pi - pi);
end

%% Modal RHA — reused from HW6 Part 4, adapted for two acceleration methods
% create vectors to store acceleration calculated for each method
% modal response = zeros(time steps, floor number, mode number)
a_modal_rel = zeros(length(ag), N, Nmodes);   % relative SDOF acc * Gamma * phi
a_modal_abs = zeros(length(ag), N, Nmodes);   % absolute SDOF acc * Gamma * phi

% load the periods, gamma, and phi from initial calculations for each mode
for n = 1:Nmodes
    T_n     = T_approx(n);
    phi_n   = PHI_approx(:,n);
    Gamma_n = Gamma_approx(n);

    % solve the sdof equation for the given mode
    % udd_abs_n = absolute (total) floor acceleration history (cm/s^2)
    [~, ~, udd_abs_n, ~,~,~,~,~] = SDOF_Response(T_n, zeta, ag, dt, 0, 0);

    % relative SDOF acceleration = absolute minus ground acceleration
    udd_rel_n = udd_abs_n - ag;

    % peak relative SDOF acc (cm/s^2)
    D_rel_max(n) = max(abs(udd_rel_n));

    % compute modal response given sdof response (history * gamma * phi)
    % gives response for mode n at all floors (y) at all times (x)
    a_modal_abs(:,:,n) = udd_abs_n * (Gamma_n * phi_n)';
    a_modal_rel(:,:,n) = udd_rel_n * (Gamma_n * phi_n)';
end
% now we have both acceleration response histories for all times,
% all floors, for all three modes

%% Sum modal contributions across 3 modes
a_floor_abs = sum(a_modal_abs, 3);   % (Nt x N) method 2 — absolute SDOF acc
a_floor_rel = sum(a_modal_rel, 3);   % (Nt x N) method 1 — relative SDOF acc

%% Recorded accelerations matrix for easy indexing
% Columns: [3rd floor, 7th floor, roof] — matches CSV columns 3-5
a_recorded = [a_3rd_rec, a_7th_rec, a_roof_rec];

%% Gradescope Q1 — peak relative SDOF accelerations in g
g = 981;   % cm/s^2
fprintf('\n--- Gradescope Q1 ---\n');
fprintf('D_rel_max Mode 1: %.4f g\n', D_rel_max(1)/g);
fprintf('D_rel_max Mode 2: %.4f g\n', D_rel_max(2)/g);
fprintf('D_rel_max Mode 3: %.4f g\n', D_rel_max(3)/g);

%% Plots — one figure per floor level
floor_idx    = [3, 7, 9];
floor_labels = {'3^{rd} Floor (j=3)', '7^{th} Floor (j=7)', 'Roof (j=9)'};

for k = 1:3
    fi = floor_idx(k);

    figure('Name', sprintf('HW7 Part A: %s', strrep(floor_labels{k}, '^{','(')));

    plot(t, a_floor_rel(:,fi), 'b-',  'LineWidth', 1.0); hold on;
    plot(t, a_floor_abs(:,fi), 'r-',  'LineWidth', 1.0);
    plot(t, a_recorded(:,k),   'k--', 'LineWidth', 1.2);

    xlabel('Time (s)');
    ylabel('Acceleration (cm/s^2)');
    title(sprintf('Floor Acceleration Time History — %s', floor_labels{k}));
    legend('Approx: Relative SDOF acc', ...
           'Approx: Absolute SDOF acc', ...
           'Recorded', ...
           'Location', 'NorthEast');
    grid on;
    xlim([t(1) t(end)]);
end

%% Save outputs for Parts B and C
save('HW7_Part1_results.mat', ...
    'T_approx', 'PHI_approx', 'Gamma_approx', ...
    'a_floor_abs', 'a_floor_rel', ...
    'a_recorded', 'ag', 'dt', 't', ...
    'N', 'H_total', 'h_story', 'x_floors', 'zeta', 'T1', 'Nmodes');