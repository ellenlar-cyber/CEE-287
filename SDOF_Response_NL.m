function [u, ud, udd_abs, Fs, Sd_inelastic, mu] = SDOF_Response_NL(Tn, z, ug, dt, u0, ud0, strength_input, input_type, method)
% Nonlinear SDOF Response - Elastoplastic (EPP) system
% Uses Newmark's method (average or linear acceleration)
% State determination based on displacement (FLAG system)
%
% INPUTS:
%   Tn             = undamped period of vibration (s)
%   z              = damping ratio
%   ug             = ground acceleration time history (cm/s^2)
%   dt             = time step of record (s)
%   u0             = initial displacement
%   ud0            = initial velocity
%   strength_input = value of either Cy or uy
%   input_type     = 'Cy' or 'uy'
%   method         = 'average' (unconditionally stable) or 'linear' (conditionally stable)
%
% OUTPUTS:
%   u             = relative displacement time history
%   ud            = relative velocity time history
%   udd_abs       = absolute acceleration time history
%   Fs            = restoring force time history
%   Sd_inelastic  = peak relative displacement (inelastic spectral ordinate)
%   mu            = peak displacement ductility demand

%% System Properties
m  = 1;
wn = 2 * pi / Tn;
k  = m * wn^2; % 2.2 Determine the tangent stiffness ki
c  = 2 * z * wn * m;
g  = 981;              % cm/s^2, must match units of ug

%% Yield Properties
% Accept either yield displacement or yield seismic coefficient
if strcmp(input_type, 'Cy')
    % Given yield seismic coefficient: Cy = Fy / (m*g)
    Fy = strength_input * m * g;
    uy = Fy / k;

elseif strcmp(input_type, 'uy')
    % Given yield displacement directly
    uy = strength_input;
    Fy = k * uy;

else
    error('input_type must be either Cy or uy')
end

%% Newmark Parameters
if strcmp(method, 'average')
    gamma = 0.5;
    beta  = 0.25;      % unconditionally stable
elseif strcmp(method, 'linear')
    gamma = 0.5;
    beta  = 1/6;       % conditionally stable, need dt < 0.551*Tn
else
    error('method must be either average or linear')
end

%% Time Step Check
% If record dt > Tn/40, interpolate the accelerogram to a finer time step
% This is needed to correctly capture yielding and unloading events
% particularly important for short period systems
dt_max = Tn / 40;
if dt > dt_max
    t_orig = (0:length(ug)-1)' * dt;
    t_fine = (0:dt_max:t_orig(end))';
    ug     = interp1(t_orig, ug(:), t_fine, 'linear');
    dt     = dt_max;
end

n = length(ug);

%% Newmark Integration Constants
% Computed once since dt is now fixed
% These appear in the effective stiffness and effective load expressions
a_const = (1/(beta*dt))*m  + (gamma/beta)*c;
b_const = (1/(2*beta))*m   + dt*(gamma/(2*beta) - 1)*c;

%% Initialize Arrays
u       = zeros(n, 1);
ud      = zeros(n, 1);
udd     = zeros(n, 1);
udd_abs = zeros(n, 1);
Fs      = zeros(n, 1);

%% Initial Conditions
u(1)   = u0;
ud(1)  = ud0;
Fs(1)  = k * u0;       % system starts elastic

% Effective earthquake force vector
p = -m * ug(:);

% Initial acceleration from equation of motion
udd(1)     = (p(1) - c*ud(1) - Fs(1)) / m;
udd_abs(1) = udd(1) + ug(1);

%% State Variables
% FLAG = 0  : elastic
% FLAG = +1 : yielding in positive direction (velocity > 0)
% FLAG = -1 : yielding in negative direction (velocity < 0)
FLAG = 0;
uyp  = +uy;   % positive yield displacement boundary
uyn  = -uy;   % negative yield displacement boundary

%% Time Stepping Loop
for i = 1:n-1
    %% Step 2.1: Effective load increment
    % Incorporates current velocity and acceleration state come back to
    % this
    dp     = p(i+1) - p(i);
    dp_hat = dp + a_const*ud(i) + b_const*udd(i);

    %% Step 2.2: Determine tangent stiffness from current FLAG
    % kt = k when elastic, kt = 0 when yielding (EPP has zero post-yield stiffness)
    if FLAG == 0
        kt = k;
    else
        kt = 0;
    end
    
    %% Step 2.3: Effective stiffness
    % Combines tangent stiffness with mass and damping contributions
    k_hat = kt + (gamma/(beta*dt))*c + (1/(beta*dt^2))*m;
    %% Step 2.4: Solve for displacement increment 
    % We understand getting this but where does the table come into play
    du = dp_hat / k_hat;

    %% Step 2.5: Update displacement and velocity
    u_new  = u(i) + du;

    ud_new = (gamma/(beta*dt))*du ...
           - (gamma/beta)*ud(i) ...
           + dt*(1 - gamma/(2*beta))*udd(i);

    %% Step 6: State Determination Based on DISPLACEMENT
    % Check displacement against yield boundaries to determine FLAG
    % Must use displacement not force so this works for bilinear systems too

    if FLAG == 0
        % Currently elastic: check if yield boundary has been crossed

        if u_new >= uyp
            % Crossed positive yield boundary
            FLAG   = +1;
            Fs_new = +Fy;

        elseif u_new <= uyn
            % Crossed negative yield boundary
            FLAG   = -1;
            Fs_new = -Fy;

        else
            % Remains elastic
            FLAG   = 0;
            Fs_new = k * u_new;
        end

    elseif FLAG == +1
        % Currently yielding positive: check for unloading via velocity reversal

        if ud_new < 0
            % Velocity has reversed: unloading event
            % Update yield boundaries: new elastic range is 2*uy wide
            % anchored at the displacement where unloading began
            uyp    = u(i);
            uyn    = u(i) - 2*uy;
            FLAG   = 0;
            Fs_new = Fy + k*(u_new - uyp);

        else
            % Continues yielding in positive direction
            Fs_new = +Fy;
        end

    else
        % FLAG == -1
        % Currently yielding negative: check for unloading via velocity reversal

        if ud_new > 0
            % Velocity has reversed: unloading event
            % Update yield boundaries
            uyn    = u(i);
            uyp    = u(i) + 2*uy;
            FLAG   = 0;
            Fs_new = -Fy + k*(u_new - uyn);

        else
            % Continues yielding in negative direction
            Fs_new = -Fy;
        end
    end

    %% Step 7: Store updated state
    u(i+1)  = u_new;
    ud(i+1) = ud(i) + ud_new; % come back to this 
    Fs(i+1) = Fs_new;

    %% Step 8: Acceleration from equation of motion
    % Must use Fs directly and NOT the incremental formula used for linear systems
    % because Fs is no longer proportional to displacement
    udd(i+1)     = (p(i+1) - c*ud(i+1) - Fs(i+1)) / m;
    udd_abs(i+1) = udd(i+1) + ug(i+1);

end

%% Output Spectral Quantities
Sd_inelastic = max(abs(u));
mu           = Sd_inelastic / uy;   % peak ductility demand

end