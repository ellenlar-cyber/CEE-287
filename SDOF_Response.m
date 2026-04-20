%% SDOF Response Function
function [u, ud, udd_abs, Sd, Sv, Sa, Spv, Spa] = SDOF_Response(Tn, z, ug, dt, u0, ud0)

% Inputs
% Tn = undamped period of vibration
% z = damping ratio
% ug = ground motion acceleration time history
% dt = time step
% u0 = initial displacement
% ud = initial velocity

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


% Define Constants of recurrance method
A = exp(-z*wn*dt) * ...
    ( (z / sqrt(1 - z^2)) * sin(wd*dt) + cos(wd*dt) );

B = exp(-z*wn*dt) * ...
    (1/wd * sin(wd*dt));

C = 1/k * ...
    (2*z / (wn*dt) + exp(-z*wn*dt)*...
    (( (1 - 2*z^2)/(wd*dt) - z/sqrt(1 - z^2) ) * sin(wd*dt) ...
        - (1 + (2*z)/(wn*dt)) * cos(wd*dt)));

D = 1/k * ...
    (1 - (2*z)/(wn*dt) + exp(-z*wn*dt) * ...
        (( (2*z^2 - 1)/(wd*dt) ) * sin(wd*dt) ...
        + ( (2*z)/(wn*dt) ) * cos(wd*dt)));

Ap = -exp(-z*wn*dt) * ...
    ( (wn / sqrt(1 - z^2)) * sin(wd*dt) );

Bp = exp(-z*wn*dt) * ...
    ( cos(wd*dt) - (z / sqrt(1 - z^2)) * sin(wd*dt) );

Cp = 1/k * ...
    (-1/(dt) + exp(-z*wn*dt) * ...
    (( wn/sqrt(1 - z^2)+(z)/(dt*sqrt(1-z^2))) * sin(wd*dt) ...
    + (1/dt) * cos(wd*dt)));

Dp = 1/(k*dt) * ...
    ( 1 - exp(-z*wn*dt)* ...
    ((z/sqrt(1-z^2)) * sin(wd*dt) + ...
    cos(wd*dt)));

% define number of iterations to perfrom
n = length(ug);

% Initialize output arrays
u = zeros(n, 1);
ud = zeros(n, 1);
udd = zeros(n,1);
udd_abs = zeros(n, 1);

% Set initial conditions
u(1) = u0;
ud(1) = ud0;

% forcing equation vector NOTE THAT IF IN TERMS OF G MUST MULTIPLY
p = -m * ug(:);

% time loop to solve for u, ud, udd
for i = 1:n-1
    % displacement calculation
    u(i+1) = A*u(i) + B*ud(i) + C*p(i) + D*p(i+1);
    % velocity calculation
    ud(i+1) = Ap*u(i) + Bp*ud(i) + Cp*p(i) + Dp*p(i+1);
    % acceleration calculation
    udd(i) = (p(i) - c*ud(i) - k*u(i))/m;
    % total acceleration calculation
    udd_abs(i) = udd(i) + ug(i);
end

% Peak relative displacement (spectral ordinate)
Sd = max(abs(u));

% Peak relative velocity (spectral ordinate)
Sv = max(abs(ud));

% Peak relative acceleration (spectral ordinate)
Sa = max(abs(udd_abs));

% Pseudo velocity spectral ordinate
Spv = wn * Sd;

% Pseduo acceleration spectral ordinate
Spa = wn^2 * Sd;

end