%% HW 2 Assignment 1 - Part A Function

function [u, ud, udd_abs, Fs, Sd, mu] = SDOF_NL(Tn, z, ug, dt, u0, ud0, uy_cy)

% Inputs
% Tn = undamped period of vibration
% z = damping ratio
% ug = ground motion acceleration time history
% dt = time step
% u0 = initial displacement
% ud = initial velocity
% uy = yield displacement
% cy = yield seismic coefficient

% u = relative displacement time history
% ud = relative velocity time history
% udd_abs = absolute acceleration time history
% Fs = restoring force time history
% Sd = inelastic disp. spectral ordinate
% mu = peak displacement ductility demand


% Define system properties
% Natural frequency
wn = 2 * pi / Tn; 
% damped frequency
wd = wn*sqrt(1-z^2);

% m 
m = 1; % assume unit mass 

% k
k = m * wn^2;

%c
c = 2 * z * wn * m;

%% Newmark Method

