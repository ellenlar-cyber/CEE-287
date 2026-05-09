%% HW6 Part 3 — Modal Response Spectrum Analysis (RSA) with SRSS
clc; clear; close all;

%% Load results from Parts 1 and 2
load('HW6_Part1_results.mat');
load('HW6_Part2_results.mat');

floors   = (1:N)';
stories  = (1:N)';
H_total  = 118 * 12;
h_story  = H_total / N;     % 157.3 in
Nmodes   = 5;               % number of modes to include

%% Step 1 — Compute modal quantities for each of the 5 modes
% Pre-allocate: rows = floors/stories, columns = modes
F_modal   = zeros(N, Nmodes);   % lateral forces  (kips)
Vs_modal  = zeros(N, Nmodes);   % story shears    (kips)
u_modal   = zeros(N, Nmodes);   % displacements   (in)
IDR_modal = zeros(N, Nmodes);   % interstory drift ratios

for n = 1:Nmodes
    phi_n   = PHI_norm(:,n);
    Gamma_n = Gamma_all(n);
    Spa_n   = Spa_modes_in(n);   % in/s^2
    Sd_n    = Sd_modes_in(n);    % in

    % Modal lateral forces at each floor
    F_modal(:,n) = m_vec .* (Gamma_n * phi_n) * Spa_n;   % kips

    % Modal story shears — sum forces from floor x to roof (BEFORE SRSS)
    for x = 1:N
        Vs_modal(x,n) = sum(F_modal(x:N, n));
    end

    % Modal floor displacements
    u_modal(:,n) = Gamma_n * phi_n * Sd_n;   % in

    % Modal interstory drift ratios
    IDR_modal(1,n)   = u_modal(1,n) / h_story;
    for x = 2:N
        IDR_modal(x,n) = (u_modal(x,n) - u_modal(x-1,n)) / h_story;
    end
end

%% Step 2 — SRSS combination (always the last step)
F_RSA   = sqrt(sum(F_modal.^2,   2));   % (Nx1)
Vs_RSA  = sqrt(sum(Vs_modal.^2,  2));   % (Nx1)
u_RSA   = sqrt(sum(u_modal.^2,   2));   % (Nx1)
IDR_RSA = sqrt(sum(IDR_modal.^2, 2));   % (Nx1)

% Base shear from SRSS of modal base shears (story 1 shear = base shear)
V_RSA  = Vs_RSA(1);
VW_RSA = V_RSA / W_tot;

%% Step 3 — Mode % contribution per floor/story
% %_n(x) = r_n(x)^2 / sum_{m=1}^{5} r_m(x)^2   (Lecture 9)

% For lateral forces
pct_F = zeros(N, 3);
sum2_F = sum(F_modal.^2, 2);   % sum of squares across 5 modes
for n = 1:3
    pct_F(:,n) = (F_modal(:,n).^2 ./ sum2_F) * 100;
end

% For story shears
pct_Vs = zeros(N, 3);
sum2_Vs = sum(Vs_modal.^2, 2);
for n = 1:3
    pct_Vs(:,n) = (Vs_modal(:,n).^2 ./ sum2_Vs) * 100;
end

% For displacements
pct_u = zeros(N, 3);
sum2_u = sum(u_modal.^2, 2);
for n = 1:3
    pct_u(:,n) = (u_modal(:,n).^2 ./ sum2_u) * 100;
end

% For IDR
pct_IDR = zeros(N, 3);
sum2_IDR = sum(IDR_modal.^2, 2);
for n = 1:3
    pct_IDR(:,n) = (IDR_modal(:,n).^2 ./ sum2_IDR) * 100;
end

%% Print results
fprintf('\n--- RSA Results (SRSS, 5 modes) ---\n');
fprintf('Base shear V_RSA   = %.2f kips\n', V_RSA);
fprintf('V_RSA / W_total    = %.4f (%.2f%%)\n', VW_RSA, VW_RSA*100);

fprintf('\n--- RSA: Floor/Story Quantities ---\n');
fprintf('%-6s %-12s %-14s %-12s %-12s\n', ...
    'Floor','F_i (kips)','Vs_x (kips)','u_i (in)','IDR_i (%%)');
for i = 1:N
    fprintf('%-6d %-12.4f %-14.4f %-12.4f %-12.4f\n', ...
        i, F_RSA(i), Vs_RSA(i), u_RSA(i), IDR_RSA(i)*100);
end

fprintf('\n--- Mode %% Contribution to Lateral Forces (modes 1-3 / 5-mode SRSS) ---\n');
fprintf('%-6s %-12s %-12s %-12s\n','Floor','Mode 1 (%%)','Mode 2 (%%)','Mode 3 (%%)');
for i = 1:N
    fprintf('%-6d %-12.1f %-12.1f %-12.1f\n', i, pct_F(i,1), pct_F(i,2), pct_F(i,3));
end

fprintf('\n--- Mode %% Contribution to Story Shears (modes 1-3 / 5-mode SRSS) ---\n');
fprintf('%-6s %-12s %-12s %-12s\n','Story','Mode 1 (%%)','Mode 2 (%%)','Mode 3 (%%)');
for i = 1:N
    fprintf('%-6d %-12.1f %-12.1f %-12.1f\n', i, pct_Vs(i,1), pct_Vs(i,2), pct_Vs(i,3));
end

fprintf('\n--- Mode %% Contribution to Displacements (modes 1-3 / 5-mode SRSS) ---\n');
fprintf('%-6s %-12s %-12s %-12s\n','Floor','Mode 1 (%%)','Mode 2 (%%)','Mode 3 (%%)');
for i = 1:N
    fprintf('%-6d %-12.1f %-12.1f %-12.1f\n', i, pct_u(i,1), pct_u(i,2), pct_u(i,3));
end

fprintf('\n--- Mode %% Contribution to IDR (modes 1-3 / 5-mode SRSS) ---\n');
fprintf('%-6s %-12s %-12s %-12s\n','Story','Mode 1 (%%)','Mode 2 (%%)','Mode 3 (%%)');
for i = 1:N
    fprintf('%-6d %-12.1f %-12.1f %-12.1f\n', i, pct_IDR(i,1), pct_IDR(i,2), pct_IDR(i,3));
end

%% Plots

% --- Figure 1: RSA profiles vs ESA ---
figure('Name','HW6 Part 3 — RSA vs ESA Profiles');

subplot(1,4,1);
barh(floors, [F_ESA, F_RSA], 0.6);
xlabel('F_i (kips)'); ylabel('Floor');
title('Lateral Forces');
legend('ESA','RSA','Location','SouthEast');
grid on; ylim([0 N+1]);

subplot(1,4,2);
plot(Vs_ESA, stories, 'b-o', 'LineWidth',1.5, 'MarkerFaceColor','b'); hold on;
plot(Vs_RSA, stories, 'r-s', 'LineWidth',1.5, 'MarkerFaceColor','r');
xlabel('V_x (kips)'); ylabel('Story');
title('Story Shears');
legend('ESA','RSA','Location','SouthEast');
grid on; ylim([0 N+1]);

subplot(1,4,3);
plot([0; u_ESA], [0; floors], 'b-o', 'LineWidth',1.5, 'MarkerFaceColor','b'); hold on;
plot([0; u_RSA], [0; floors], 'r-s', 'LineWidth',1.5, 'MarkerFaceColor','r');
xlabel('u_i (in)'); ylabel('Floor');
title('Displacements');
legend('ESA','RSA','Location','SouthEast');
grid on; ylim([0 N+1]);

subplot(1,4,4);
plot(IDR_ESA*100, stories, 'b-o', 'LineWidth',1.5, 'MarkerFaceColor','b'); hold on;
plot(IDR_RSA*100, stories, 'r-s', 'LineWidth',1.5, 'MarkerFaceColor','r');
xlabel('IDR (%)'); ylabel('Story');
title('Interstory Drift Ratios');
legend('ESA','RSA','Location','SouthEast');
grid on; ylim([0 N+1]);

sgtitle('RSA (5 modes, SRSS) vs ESA');

% --- Figure 2: Mode % contributions ---
figure('Name','HW6 Part 3 — Mode Contributions');

subplot(1,4,1);
plot(pct_F(:,1), floors,'b-o', pct_F(:,2), floors,'r-s', pct_F(:,3), floors,'g-^','LineWidth',1.5);
xlabel('Contribution (%)'); ylabel('Floor');
title('Forces'); legend('Mode 1','Mode 2','Mode 3'); grid on; ylim([0 N+1]);

subplot(1,4,2);
plot(pct_Vs(:,1),stories,'b-o', pct_Vs(:,2),stories,'r-s', pct_Vs(:,3),stories,'g-^','LineWidth',1.5);
xlabel('Contribution (%)'); ylabel('Story');
title('Story Shears'); legend('Mode 1','Mode 2','Mode 3'); grid on; ylim([0 N+1]);

subplot(1,4,3);
plot(pct_u(:,1),floors,'b-o', pct_u(:,2),floors,'r-s', pct_u(:,3),floors,'g-^','LineWidth',1.5);
xlabel('Contribution (%)'); ylabel('Floor');
title('Displacements'); legend('Mode 1','Mode 2','Mode 3'); grid on; ylim([0 N+1]);

subplot(1,4,4);
plot(pct_IDR(:,1),stories,'b-o', pct_IDR(:,2),stories,'r-s', pct_IDR(:,3),stories,'g-^','LineWidth',1.5);
xlabel('Contribution (%)'); ylabel('Story');
title('IDR'); legend('Mode 1','Mode 2','Mode 3'); grid on; ylim([0 N+1]);

sgtitle('Mode % Contributions to Response Quantities (Modes 1-3 vs 5-mode SRSS)');

%% Save for Part 7
save('HW6_Part3_results.mat', ...
     'F_RSA','Vs_RSA','u_RSA','IDR_RSA','V_RSA','VW_RSA', ...
     'F_modal','Vs_modal','u_modal','IDR_modal', ...
     'pct_F','pct_Vs','pct_u','pct_IDR');