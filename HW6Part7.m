%% HW 6 Part 7
clc; clear; close all;

%% Load results from Parts 1,2,3,4,5,6
load('HW6_Part1_results.mat');
load('HW6_Part2_results.mat');
load('HW6_Part3_results.mat');
load('HW6_Part4_results.mat');
load('HW6_Part5_results.mat');

floors  = (1:N)';
stories = (1:N)';

%  SUMMARY TABLES

fprintf('PART 7: Summary tables\n');

% displacements
fprintf('\nFloor Displacements u_i (in)\n');
fprintf('%-6s %-10s %-10s %-10s %-10s %-10s\n', ...
    'Floor','ESA','RSA','RHA','Approx','Exact');
fprintf('%s\n', repmat('-',1,56));
for i = 1:N
    if i == N   % only print approx/exact at roof (floor 9)
        fprintf('%-6d %-10.4f %-10.4f %-10.4f %-10.4f %-10.4f\n', ...
            i, u_ESA(i), u_RSA(i), u_RHA(i), u_roof_approx, u_roof_exact);
    else
        fprintf('%-6d %-10.4f %-10.4f %-10.4f %-10s %-10s\n', ...
            i, u_ESA(i), u_RSA(i), u_RHA(i), '---', '---');
    end
end

% idr
fprintf('\nInterstory Drift Ratios IDR_i (%%)\n');
fprintf('%-6s %-10s %-10s %-10s %-10s %-10s\n', ...
    'Story','ESA','RSA','RHA','Approx','Exact');
fprintf('%s\n', repmat('-',1,56));

% find which story has peak IDR for RHA (to know where to print approx)
[~, peak_IDR_story] = max(IDR_RHA);

for i = 1:N
    if i == peak_IDR_story
        fprintf('%-6d %-10.4f %-10.4f %-10.4f %-10.4f %-10.4f\n', ...
            i, IDR_ESA(i)*100, IDR_RSA(i)*100, IDR_RHA(i)*100, ...
            IDR_approx*100, IDR_exact*100);
    else
        fprintf('%-6d %-10.4f %-10.4f %-10.4f %-10s %-10s\n', ...
            i, IDR_ESA(i)*100, IDR_RSA(i)*100, IDR_RHA(i)*100, '---', '---');
    end
end

% lateral forces
fprintf('\nLateral Forces F_i (kips)\n');
fprintf('%-6s %-10s %-10s %-10s\n', 'Floor','ESA','RSA','RHA');
fprintf('%s\n', repmat('-',1,36));
for i = 1:N
    fprintf('%-6d %-10.4f %-10.4f %-10.4f\n', ...
        i, F_ESA(i), F_RSA(i), F_RHA(i));
end

% story shears
fprintf('\nStory Shears Vs_x (kips)\n');
fprintf('%-6s %-10s %-10s %-10s\n', 'Story','ESA','RSA','RHA');
fprintf('%s\n', repmat('-',1,36));
for i = 1:N
    fprintf('%-6d %-10.4f %-10.4f %-10.4f\n', ...
        i, Vs_ESA(i), Vs_RSA(i), Vs_RHA(i));
end


% base shear
fprintf('\nBase Shear Summary\n');
fprintf('%-10s %-12s %-12s\n','Method','V (kips)','V/W (%%)');
fprintf('%s\n', repmat('-',1,36));
fprintf('%-10s %-12.2f %-12.4f\n', 'ESA', V_ESA, VW_ESA*100);
fprintf('%-10s %-12.2f %-12.4f\n', 'RSA', V_RSA, VW_RSA*100);
fprintf('%-10s %-12.2f %-12.4f\n', 'RHA', V_RHA, VW_RHA*100);

%% comparison plots for ESA, RSA, RHA
figure('Name','HW6 Part 7: All Response Quantities');

% displacements
subplot(2,2,1);
plot([0;u_ESA], [0;floors], 'b-o', 'LineWidth',1.5); hold on;
plot([0;u_RSA], [0;floors], 'r-o', 'LineWidth',1.5);
plot([0;u_RHA], [0;floors], 'g-o', 'LineWidth',1.5);
xlabel('u_i (in)'); ylabel('Floor');
title('Displacements');
legend('ESA','RSA','RHA','Location','SouthEast');
grid on; ylim([0 N+1]);

% idr
subplot(2,2,2);
plot(IDR_ESA*100, stories, 'b-o', 'LineWidth',1.5); hold on;
plot(IDR_RSA*100, stories, 'r-o', 'LineWidth',1.5);
plot(IDR_RHA*100, stories, 'g-o', 'LineWidth',1.5);
xlabel('IDR (%)'); ylabel('Story');
title('Interstory Drift Ratios');
legend('ESA','RSA','RHA','Location','NorthEast');
grid on; ylim([0 N+1]);

% lateral forces
subplot(2,2,4);
plot(F_ESA, floors, 'b-o', 'LineWidth',1.5); hold on;
plot(F_RSA, floors, 'r-o', 'LineWidth',1.5);
plot(F_RHA, floors, 'g-o', 'LineWidth',1.5);
xlabel('F_i (kips)'); ylabel('Floor');
title('Lateral Forces');
legend('ESA','RSA','RHA','Location','SouthEast');
grid on; ylim([0 N+1]);

% story shear
subplot(2,2,3);
plot(Vs_ESA, stories, 'b-o', 'LineWidth',1.5); hold on;
plot(Vs_RSA, stories, 'r-o', 'LineWidth',1.5);
plot(Vs_RHA, stories, 'g-o', 'LineWidth',1.5);
xlabel('V_x (kips)'); ylabel('Story');
title('Story Shears');
legend('ESA','RSA','RHA','Location','SouthEast');
grid on; ylim([0 N+1]);

sgtitle('Part 7: ESA vs RSA vs RHA: All Response Quantities');

%% comparison of mode 1 contribution to base shear
% Mode 1 base shear from RSA (already in Part 3)
V_mode1_RSA = Vs_modal(1,1);
pct_mode1_RSA = pct_Vs(1,1);

% Mode 1 base shear from RHA (already in Part 4)
V_mode1_RHA = abs(Vs_modes_atPeak(1,1));
pct_mode1_RHA = pct_Vs_RHA(1,1);


fprintf('\nMode 1 Contribution to Base Shear\n');
fprintf('ESA:  V = %.2f kips (100%% by definition)\n', V_ESA);
fprintf('RSA:  Mode 1 = %.2f kips = %.1f%% of SRSS base shear\n', ...
    V_mode1_RSA, pct_mode1_RSA);
fprintf('RHA:  Mode 1 = %.2f kips = %.1f%% at peak base shear instant\n', ...
    V_mode1_RHA, pct_mode1_RHA);