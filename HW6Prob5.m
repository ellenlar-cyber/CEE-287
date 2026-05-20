%% HW 6 part 5
clc; clear; close all;

%% Load results from Parts 1 and 2
load('HW6_Part1_results.mat');
load('HW6_Part2_results.mat');

floors   = (1:N)';
stories  = (1:N)';
H_total  = 118 * 12;
h_story  = H_total / N;     % 157.3 in
Nmodes   = 5;               % number of modes to include

% pull max sd for 1st mode from earlier parts (part 2)
Sd1     = Sd_modes_in(1);

n = 1;
% calculation of gamma for each mode
gamma1 = 4 * (-1) ^ (n-1) / (2*n*pi - pi);
% Define the height vector for mode shape calculation
x = linspace(h_story, H_total, N);
% calculate phi for each story
phi1(:) = (-1)^(n-1) * sin( ((2*n -1) * pi * x) / (2*H_total));

% calculate disp for each story for each mode
% Calculate displacement for each story for each mode
u_story(:) = Sd1 * (gamma1 * phi1);
%% Approximate calculation
% beta1 approx is gamma 1 times phi 1
beta1_approx = gamma1;
% u roof is disp at story 9 (top)
u_roof_approx = (u_story(9));
% beta 2 from lecture slides
beta2_approx = 1.505; % from lecture slides
% calculate interstory drift ratio
IDR_approx = beta2_approx*u_roof_approx/(H_total);

%% "Exact" Calculation
% for non approximate calculation, use beta 1 from gamma and phi
% calcualtions
% beta 1 will be effective mode shape 1
beta1_exact = Beta1_1;
%uroof
u_roof_exact = beta1_exact* Sd1;
% pull Beta 2 from part 2
beta2_exact = Beta2_1;
% interstory drift ratio
IDR_exact = Sd1*beta1_exact*beta2_exact/(H_total);

fprintf('     PART 5: APPROXIMATE METHOD RESULTS    \n');
fprintf('Sd(T1) = %.4f in\n\n', Sd1);

fprintf('%-22s %12s %12s\n', 'Parameter', 'Approx', 'Exact');
fprintf('%-22s %12.4f %12.4f\n', 'beta1',        beta1_approx,   beta1_exact);
fprintf('%-22s %12.4f %12.4f\n', 'beta2 ', beta2_approx,   beta2_exact);
fprintf('%-22s %12.4f %12.4f\n', 'u_roof (in)',  u_roof_approx,  u_roof_exact);
fprintf('%-22s %12.4f %12.4f\n', 'Peak IDR (%)', IDR_approx*100, IDR_exact*100);

save('HW6_Part5_results.mat', ...
     'beta1_approx','beta1_exact','beta2_approx','beta2_exact', ...
     'u_roof_approx','u_roof_exact','IDR_approx','IDR_exact')