%% Hw 6 part 4
clc; clear; close all;

%% Load results from Parts 1 and 2
load('HW6_Part1_results.mat');
load('HW6_Part2_results.mat');

floors   = (1:N)';
stories  = (1:N)';
H_total  = 118 * 12;
h_story  = H_total / N;     % 157.3 in
Nmodes   = 5;               % number of modes to include

% number of computations to perform is length of ground acceleration
Ncalc = length(ag);

% create vectors to store displacement and accelertation calculated
% modal response = zeros(time of calculation, floor number, mode number)
u_modalResponse = zeros(length(ag), length(stories), Nmodes);
a_modalResponse = zeros(length(ag), length(stories), Nmodes);

% load the periods, gamma, and phi from initial calculations for each mode
for n = 1:Nmodes
    T_n  = R.T(n);
    phi_n   = PHI_norm(:,n);
    Gamma_n = Gamma_all(n);

    % solve the sdof equation for the given mode
    [u_hist,~,udd_hist, ~,~,~,~, ~] = SDOF_Response(T_n, zeta, ag, dt, 0, 0);

    % convert from cm to in
    u_hist = u_hist/2.54;
    udd_hist = udd_hist/2.54;

    % compute modal response given sdof response (history*gamma*phi)
    % will give reponse for given node n at all floors (y) at all times (x)
    u_modalResponse(:, :, n) = u_hist * (Gamma_n * phi_n)'; 
    a_modalResponse(:, :, n) = udd_hist * (Gamma_n*phi_n)';
end
% now we have the displacement and acceleration response history for all
% times, all floors, for all modes

%% Calculate peak for each mode independently

% Compute summed reponse including all modes
% initialize story shears and interstory drift arrays for calculation
% across all time steps
Vs_hist = zeros(length(ag), N);
IDR_hist = zeros(length(ag), N);

% sum modal contributions across 5 modes
u_hist_total = sum(u_modalResponse(:,:,1:5), 3);   % (Ntimes x N)
a_hist_total = sum(a_modalResponse(:,:,1:5), 3);   % (Ntimes x N)

% Lateral forces at each floor: f = a *m
F_hist = a_hist_total .* m_vec';  

% Modal story shears — sum forces from floor x to roof (BEFORE SRSS)
for x = 1:N
    Vs_hist(:,x) = sum(F_hist(:, x:N), 2);   
end
    
% Interstory drift ratios
IDR_hist(:,1) = u_hist_total(:,1) / h_story;  
for x = 2:N
    IDR_hist(:,x) = (u_hist_total(:,x) - u_hist_total(:,x-1)) / h_story;
end

% peak is max absolute value over entire reposne history analysis
u_RHA   = max(abs(u_hist_total),[], 1)';
Vs_RHA  = max(abs(Vs_hist),[], 1)';
F_RHA   = max(abs(F_hist),[], 1)';
IDR_RHA = max(abs(IDR_hist),[], 1)';

% base shear will be the first entry of story shear (at floor 1)
V_RHA = Vs_RHA(1);

% as fraction of total weight of building:
VW_RHA = V_RHA / W_tot;

fprintf('\n--- RHA Results (5 modes) ---\n');
fprintf('Base shear V_RHA   = %.2f kips\n', V_RHA);
fprintf('V_RHA / W_total    = %.4f (%.2f%%)\n', VW_RHA, VW_RHA*100);

fprintf('\n%-6s %-12s %-14s %-12s %-12s\n', ...
    'Floor','F_i (kips)','Vs_x (kips)','u_i (in)','IDR_i (%%)');
for i = 1:N
    fprintf('%-6d %-12.4f %-14.4f %-12.4f %-12.4f\n', ...
        i, F_RHA(i), Vs_RHA(i), u_RHA(i), IDR_RHA(i)*100);
end

%% compute the index where peak values occur
idx_u   = zeros(N, 1);   % time index of peak displacement at each floor
idx_F   = zeros(N, 1);   % time index of peak force at each floor
idx_Vs  = zeros(N, 1);   % time index of peak shear at each story
idx_IDR = zeros(N, 1);   % time index of peak IDR at each story

for i = 1:N
    idx_u(i)= find(abs(u_hist_total(:,i)) == u_RHA(i),1);
    idx_F(i)= find(abs(F_hist(:,i))== F_RHA(i),1);
    idx_Vs(i)= find(abs(Vs_hist(:,i))== Vs_RHA(i),1);
    idx_IDR(i)= find(abs(IDR_hist(:,i))== IDR_RHA(i),1);
end

%% Extract modal contributions at each floor's own peak time
% Result matrices: rows = floors, columns = modes
u_modes_atPeak   = zeros(N, Nmodes);
F_modes_atPeak   = zeros(N, Nmodes);
Vs_modes_atPeak  = zeros(N, Nmodes);
IDR_modes_atPeak = zeros(N, Nmodes);

for i = 1:N
    for n = 1:Nmodes
        % displacement of floor i from mode n, at the time floor i peaks
        u_modes_atPeak(i,n) = u_modalResponse(idx_u(i), i, n);

        % force at floor i from mode n, at the time floor i's force peaks
        F_modes_atPeak(i,n) = a_modalResponse(idx_F(i), i, n) * m_vec(i);

        % story shear at story i from mode n, at the time story i peaks
        
        F_at_t = a_modalResponse(idx_Vs(i), :, n)' .* m_vec; 
        % shear = sum of forces from floor i to roof at that instant
        Vs_modes_atPeak(i,n) = sum(F_at_t(i:N));

        % IDR at story i from mode n, at the time story i's IDR peaks
        if i == 1
            IDR_modes_atPeak(i,n) = u_modalResponse(idx_IDR(i), 1, n) / h_story;
        else
            IDR_modes_atPeak(i,n) = (u_modalResponse(idx_IDR(i), i,   n) - ...
                                     u_modalResponse(idx_IDR(i), i-1, n)) / h_story;
        end
    end
end

%% Modal percentage contributions at each floor's own peak
pct_u = (abs(u_modes_atPeak)./ sum(abs(u_modes_atPeak),2)) * 100;
pct_F = (abs(F_modes_atPeak)./ sum(abs(F_modes_atPeak),2)) * 100;
pct_Vs_RHA = (abs(Vs_modes_atPeak)./ sum(abs(Vs_modes_atPeak),2)) * 100;
pct_IDR = (abs(IDR_modes_atPeak)./ sum(abs(IDR_modes_atPeak),2)) * 100;

%% Print modal contributions at peak
fprintf('\n========================================\n');
fprintf('   MODAL CONTRIBUTIONS AT PEAK VALUES   \n');
fprintf('========================================\n');

fprintf('\n--- Peak Displacements (in) ---\n');
fprintf('%-6s %-8s %-8s %-8s %-8s %-8s %-8s\n', ...
    'Floor','M1','M2','M3','M4','M5','Total');
for i = 1:N
    fprintf('%-6d %-8.4f %-8.4f %-8.4f %-8.4f %-8.4f %-8.4f\n', ...
        i, u_modes_atPeak(i,:), u_RHA(i));
end

fprintf('\n--- Peak Lateral Forces (kips) ---\n');
fprintf('%-6s %-8s %-8s %-8s %-8s %-8s %-8s\n', ...
    'Floor','M1','M2','M3','M4','M5','Total');
for i = 1:N
    fprintf('%-6d %-8.4f %-8.4f %-8.4f %-8.4f %-8.4f %-8.4f\n', ...
        i, F_modes_atPeak(i,:), F_RHA(i));
end

fprintf('\n--- Peak Story Shears (kips) ---\n');
fprintf('%-6s %-8s %-8s %-8s %-8s %-8s %-8s\n', ...
    'Story','M1','M2','M3','M4','M5','Total');
for i = 1:N
    fprintf('%-6d %-8.4f %-8.4f %-8.4f %-8.4f %-8.4f %-8.4f\n', ...
        i, Vs_modes_atPeak(i,:), Vs_RHA(i));
end

fprintf('\n--- Peak IDR (%%) ---\n');
fprintf('%-6s %-8s %-8s %-8s %-8s %-8s %-8s\n', ...
    'Story','M1','M2','M3','M4','M5','Total');
for i = 1:N
    fprintf('%-6d %-8.4f %-8.4f %-8.4f %-8.4f %-8.4f %-8.4f\n', ...
        i, IDR_modes_atPeak(i,:)*100, IDR_RHA(i)*100);
end

fprintf('\n--- Mode %% Contributions at Peak ---\n');
fprintf('\nDisplacements:\n');
fprintf('%-6s %-8s %-8s %-8s %-8s %-8s\n','Floor','M1%%','M2%%','M3%%','M4%%','M5%%');
for i = 1:N
    fprintf('%-6d %-8.1f %-8.1f %-8.1f %-8.1f %-8.1f\n', i, pct_u(i,:));
end

fprintf('\nStory Shears:\n');
fprintf('%-6s %-8s %-8s %-8s %-8s %-8s\n','Story','M1%%','M2%%','M3%%','M4%%','M5%%');
for i = 1:N
    fprintf('%-6d %-8.1f %-8.1f %-8.1f %-8.1f %-8.1f\n', i, pct_Vs_RHA(i,:));
end


%% HW 6 Problem 6
% (a) plot lateral forces when roof displacement is maximum
% (b) plot lateral forces when base shear is maximum
% (a) plot shear forces when roof displacement is maximum
% (b) plot shear forces when base shear is maximum

% from part 4, response history analysis
% find max roof disp
% find max base shear
% these were done in part 4 using combined response history analysis

% pull index of peak roof disp value - will be i = floor 9
idx_u= find(abs(u_hist_total(:,9)) == u_RHA(9),1);
% pull index of peak base shear value - will be i = floor 1
idx_Vs= find(abs(Vs_hist(:,1))== Vs_RHA(1),1);

% compute lateral forces when roof disp is max
% compute at time associated with index, for all floors
F_roofmax  = abs(F_hist(idx_u,  :))';       

% compute lateral forces when base shear is max
% compute at time associated with index, for all floors
F_basemax  = abs(F_hist(idx_Vs,  :))';

% get lateral forces from ESA (part 2)
F_esa = F_ESA;
% get lateral forces maximum calculated for RHA (part 4)
F_max = F_RHA;

% SIMILAR BUT NOW SHEAR FORCES
% for calcualtion of shear force when roof disp is max
Vs_roofmax = abs(Vs_hist(idx_u, :))';
% for calcualtion of shear force when base shear is max
Vs_basemax = abs(Vs_hist(idx_Vs, :))';

% get shear forces from ESA (part 2)
Vs_esa = Vs_ESA;
% get shear forces maximum calculated for RHA (part 4)
Vs_max = Vs_RHA;


%% Figure — Part 6
figure('Name', 'Part 6');

% --- Lateral Forces ---
subplot(1,2,1);
plot(F_esa,     floors, 'b-o', 'LineWidth', 1.5); hold on;
plot(F_roofmax, floors, 'r-o', 'LineWidth', 1.5);
plot(F_basemax, floors, 'g-o', 'LineWidth', 1.5);
plot(F_max,     floors, 'k-o', 'LineWidth', 1.5);
xlabel('Lateral Force (kips)');
ylabel('Floor');
title('Lateral Forces');
legend('ESA', 'Peak Roof Disp', 'Peak Base Shear', 'RHA Peak', ...
    'Location', 'SouthEast');
grid on; ylim([0 N+1]);

% --- Story Shears ---
subplot(1,2,2);
plot(Vs_esa,     stories, 'b-o', 'LineWidth', 1.5); hold on;
plot(Vs_roofmax, stories, 'r-o', 'LineWidth', 1.5);
plot(Vs_basemax, stories, 'g-o', 'LineWidth', 1.5);
plot(Vs_max,     stories, 'k-o', 'LineWidth', 1.5);
xlabel('Story Shear (kips)');
ylabel('Story');
title('Story Shears');
legend('ESA', 'Peak Roof Disp', 'Peak Base Shear', 'RHA Peak', ...
    'Location', 'SouthEast');
grid on; ylim([0 N+1]);

sgtitle('Part 6: Force Distributions');

save('HW6_Part4_results.mat', ...
     'u_RHA', 'F_RHA', 'Vs_RHA', 'IDR_RHA', ...
     'V_RHA', 'VW_RHA', ...
     'u_hist_total', 'F_hist', 'Vs_hist', 'IDR_hist', ...
     'u_modes_atPeak', 'F_modes_atPeak', 'Vs_modes_atPeak', 'IDR_modes_atPeak', ...
     'pct_u', 'pct_F', 'pct_Vs_RHA', 'pct_IDR', ...
     'F_roofmax', 'F_basemax', 'Vs_roofmax', 'Vs_basemax');