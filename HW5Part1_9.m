%% CEE-287 Assignment 5 — MDOF Sensitivity Analysis
%  9-story shear building eigenvalue analysis and sensitivity studies.
%  Requires: mdof_analysis.m, stiffness_vector.m, mass_vector.m

clc; clear; close all;

%  PROBLEM 1 — Baseline: 9-story, uniform m and k
% Use unit values; period ratios and mode shapes are scale-independent.
N = 9;
m = 1;   % floor mass
k = 1;   % story stiffness

m_vec = m * ones(N, 1);
k_vec = k * ones(N, 1);

fprintf('\nPROBLEM 1: Baseline (N=9, uniform m, uniform k)\n');
R1 = mdof_analysis(N, m_vec, k_vec, 'Problem 1 — Baseline');

%  PROBLEM 2 — Double both m and k (2m, 2k)
fprintf('\nPROBLEM 2: 2m, 2k\n');
R2 = mdof_analysis(N, 2*m_vec, 2*k_vec, 'Problem 2 — 2m, 2k');

%  PROBLEM 3 — Double mass only (2m, k)
% T1 increases by sqrt(2); period RATIOS and mode shapes are unchanged.
fprintf('\nPROBLEM 3: 2m, k\n');
R3 = mdof_analysis(N, 2*m_vec, k_vec, 'Problem 3 — 2m, k');


%  PROBLEMS 4-6 — Variable stiffness: k(j) = 1-(1-delta)*(j/N)^lambda

% Problem 4: delta=0.75, lambda=1  (mild linear reduction)
fprintf('\nPROBLEM 4: Variable k, delta=0.75, lambda=1\n');
k4 = stiffness_vector(N, 0.75, 1);
R4 = mdof_analysis(N, m_vec, k4, 'Problem 4 — delta=0.75, lam=1');

% Problem 5: delta=0.2, lambda=1  (steep linear reduction)
fprintf('\nPROBLEM 5: Variable k, delta=0.2, lambda=1\n');
k5 = stiffness_vector(N, 0.2, 1);
R5 = mdof_analysis(N, m_vec, k5, 'Problem 5 — delta=0.2, lam=1');

% Problem 6: delta=0.2, lambda=2  (parabolic reduction)
fprintf('\nPROBLEM 6: Variable k, delta=0.2, lambda=2\n');
k6 = stiffness_vector(N, 0.2, 2);
R6 = mdof_analysis(N, m_vec, k6, 'Problem 6 — delta=0.2, lam=2');

%  PROBLEM 7 — Parabolic mass reduction, uniform k (delta=0.2, lambda=2)
fprintf('\nPROBLEM 7: Variable m (delta=0.2, lam=2), uniform k\n');
m7 = mass_vector(N, 0.2, 2);
R7 = mdof_analysis(N, m7, k_vec, 'Problem 7 — variable m, uniform k');

%  PROBLEM 8 — 20-story building, parabolic k (delta=0.2, lambda=2)
fprintf('\nPROBLEM 8: N=20, delta=0.2, lambda=2\n');
N8 = 20;
m8 = ones(N8, 1);
k8 = stiffness_vector(N8, 0.2, 2);
R8 = mdof_analysis(N8, m8, k8, 'Problem 8 — N=20, delta=0.2, lam=2');

%  PROBLEM 9 — Soft story at ground level (9-story)
fprintf('\nPROBLEM 9: Soft story k(1)=0.63\n');
k9 = stiffness_vector(N, 0.2, 2);
k9(1) = 0.63;   % override ground story to create soft story irregularity
R9 = mdof_analysis(N, m_vec, k9, 'Problem 9 — Soft story k(1)=0.63');

%  SUMMARY TABLES
fprintf('\n Period Ratio Summary Table \n');
fprintf('%-22s  T2/T1    T3/T1\n', 'Problem');
fprintf('%-22s  %.4f   %.4f\n', 'P1 (uniform)',    R1.T2T1, R1.T3T1);
fprintf('%-22s  %.4f   %.4f\n', 'P4 (d=0.75,l=1)', R4.T2T1, R4.T3T1);
fprintf('%-22s  %.4f   %.4f\n', 'P5 (d=0.2,l=1)',  R5.T2T1, R5.T3T1);
fprintf('%-22s  %.4f   %.4f\n', 'P6 (d=0.2,l=2)',  R6.T2T1, R6.T3T1);
fprintf('%-22s  %.4f   %.4f\n', 'P7 (var mass)',    R7.T2T1, R7.T3T1);
fprintf('%-22s  %.4f   %.4f\n', 'P8 (N=20)',        R8.T2T1, R8.T3T1);
fprintf('%-22s  %.4f   %.4f\n', 'P9 (soft story)',  R9.T2T1, R9.T3T1);
