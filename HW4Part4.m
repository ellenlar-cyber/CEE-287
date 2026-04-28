%% HW 4 Number 4
clc; clear; close all

%% BILINEAR FUNCTION WITH AUTOMATIC TIME STEPPING COMMENTED OUT
function [u, ud, udd_abs, Fs, Sd_inelastic, mu] = Bilinear_SDOF_Response_NL(Tn, z, ug, dt, u0, ud0, strength_input, input_type,alpha, method)
% Bilinear SDOF system
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
%   alpha          = post-yield stiffness ratio 
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

% Add initial biliateral inital stifness and strength
k_epp  = k * (1 - alpha);    % EPP component initial stiffness
Fy_epp = Fy * (1 - alpha);   % EPP component yield force
% Note: uy = Fy/k is unchanged — the yield displacement of the full system is the same

%% Newmark Parameters
if strcmp(method, 'average')
    gamma = 1/2;
    beta  = 1/4;      % unconditionally stable
elseif strcmp(method, 'linear')
    gamma = 1/2;
    beta  = 1/6;       % conditionally stable, need dt < 0.551*Tn
else
    error('method must be either average or linear')
end

%% Time Step Check
% If record dt > Tn/40, interpolate the accelerogram to a finer time step
% This is needed to correctly capture yielding and unloading events
% particularly important for short period systems
% dt_max = Tn / 40;
% if dt > dt_max
%     t_orig = (0:length(ug)-1)' * dt;
%     t_fine = (0:dt_max:t_orig(end))';
%     ug     = interp1(t_orig, ug(:), t_fine, 'linear');
%     dt     = dt_max;
% end

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
        kt = alpha * k; % bilinear: post-yield stiffness from linear spring component
    end
    
    %% Step 2.3: Effective stiffness
    % Combines tangent stiffness with mass and damping contributions
    k_hat = kt + (gamma/(beta*dt))*c + (1/(beta*dt^2))*m;
    %% Step 2.4: Solve for displacement increment 
    % We understand getting this but where does the table come into play
    du = dp_hat / k_hat;

    %% Step 2.5: Update displacement and velocity
    u_new  = u(i) + du; % eq 2.7

    ud_inc = (gamma/(beta*dt))*du ...
           - (gamma/beta)*ud(i) ...
           + dt*(1 - gamma/(2*beta))*udd(i);
    ud_new = ud_inc + ud(i); % eq 2.7

    %% Step 6: State Determination Based on DISPLACEMENT
    % Check displacement against yield boundaries to determine FLAG
    % Must use displacement not force so this works for bilinear systems too

    if FLAG == 0
        % Currently elastic: check if yield boundary has been crossed

        if u_new >= uyp
            % Crossed positive yield boundary
            FLAG   = +1; % when the Flag is changed it goes to thier repsective elseif statment to get correct Fs_new
            Fs_new = +Fy_epp + (alpha *k * u_new);

        elseif u_new <= uyn
            % Crossed negative yield boundary
            FLAG   = -1; % when the Flag is changed it goes to thier repsective elseif statment to get correct Fs_new
            Fs_new = -Fy_epp + (alpha *k * u_new);

        else
            % Remains elastic
            FLAG   = 0;
            Fs_new =  Fs(i) + k*(u_new - u(i));
        end

    elseif FLAG == +1
        % Currently yielding positive: check for unloading via velocity reversal

        if ud_new < 0
            % Velocity has reversed: unloading event
            % Update yield boundaries: new elastic range is 2*uy wide
            % anchored at the displacement where unloading began
            uyp    = u(i);
            uyn    = u(i) - 2*uy;
            FLAG   = 0; % needs to change flag back to 0 so when it starts the next iteration it goes to the Flag ==0 if statement that will then send to correct Fs calc depending on the new time step conditions
            Fs_new = Fy_epp + alpha*k*u(i);

        else
            % Continues yielding in positive direction
            Fs_new = +Fy_epp + alpha*k*u_new;
        end

    else
        % FLAG == -1
        % Currently yielding negative: check for unloading via velocity reversal

        if ud_new > 0
            % Velocity has reversed: unloading event
            % Update yield boundaries
            uyn    = u(i);
            uyp    = u(i) + 2*uy;
            FLAG   = 0; % needs to change flag back to 0 so when it starts the next iteration it goes to the Flag ==0 if statement that will then send to correct Fs calc depending on the new time step conditions
            Fs_new = -Fy_epp + alpha*k*u(i);

        else
            % Continues yielding in negative direction
            Fs_new = -Fy_epp + alpha*k*u_new;
        end
    end

    %% Step 7: Store updated state
    u(i+1)  = u_new;
    ud(i+1) = ud_new; % come back to this 
    Fs(i+1) = Fs_new;

    %% Step 8: Acceleration from equation of motion
    % Must use Fs directly and NOT the incremental formula used for linear systems
    % because Fs is no longer proportional to displacement
    % Recommended equation to use not the one listed 2.7 for relative
    % acceleration
    udd(i+1)     = (p(i+1) - c*ud(i+1) - Fs(i+1)) / m;
    udd_abs(i+1) = udd(i+1) + ug(i+1);

    % Reset FLAG to 0
    %FLAG = 0;

end

%% Output Spectral Quantities
% Inelastic displacement spectral ordinate calculation (peak relative
        % displacement)
Sd_inelastic = max(abs(u));

% Peak displacement ductility demand (peak relative displacement/yield
        % displacement)
mu           = Sd_inelastic / uy;   % peak ductility demand

end



%% Define structure conditions
T = 0.4;
z = 0.02;
m = 1;
Cy = 0.4;
a = 0.0;

%% Load El Centro Record
% Data has two columns: time (s) and acceleration (g)
data = readmatrix('El_Centro_NS.xls');
t = data(:, 1);
ag_g     = data(:, 2);       % acceleration in g

% Convert 
g   = 386.2;                   
ag  = ag_g * g;              

%% For time step 0.02 case (record)
dt  = t(2) - t(1);   

% Call function
[u, ~, ~, ~, Sd, ~] = Bilinear_SDOF_Response_NL(T, z, ag, dt, 0, 0, Cy, 'Cy',a, 'linear');

% Store values
t_02 = t;
u_02 = u;
Sd_02 = Sd;
u_res_02 = u(end);  % check this for residual disp calc

%% For time step case 0.0100 

dt = 0.0100;
t_new = (t(1): dt : t(end));
ag_new = interp1(t, ag, t_new, 'linear');

% Call function
[u, ~, ~, ~, Sd, ~] = Bilinear_SDOF_Response_NL(T, z, ag_new, dt, 0, 0, Cy, 'Cy',a, 'linear');

% Store values
t_01 = t_new;
u_01 = u;
Sd_01 = Sd;
u_res_01 = u(end);  % check this for residual disp calc

%% For time step case 0.0020

dt = 0.0020;
t_new = (t(1): dt : t(end));
ag_new = interp1(t, ag, t_new, 'linear');

% Call function
[u, ~, ~, ~, Sd, ~] = Bilinear_SDOF_Response_NL(T, z, ag_new, dt, 0, 0, Cy, 'Cy',a, 'linear');

% Store values
t_002 = t_new;
u_002 = u;
Sd_002 = Sd;
u_res_002 = u(end);  % check this for residual disp calc

%% For time step case 0.0005

dt = 0.0005;
t_new = (t(1): dt : t(end));
ag_new = interp1(t, ag, t_new, 'linear');

% Call function
[u, ~, ~, ~, Sd, ~] = Bilinear_SDOF_Response_NL(T, z, ag_new, dt, 0, 0, Cy, 'Cy',a, 'linear');

% Store values
t_0005 = t_new;
u_0005 = u;
Sd_0005 = Sd;
u_res_0005 = u(end);  % check this for residual disp calc

%% Print Results
%% Print Results Table

fprintf('Problem 4: Bilinear SDOF El Centro NS\n');
fprintf('%-12s %-18s %-18s\n', 'dt (s)', 'Peak Sd (in)', 'Residual u (in)');
fprintf('%-12.4f %-18.4f %-18.4f\n', 0.0200, Sd_02,   u_res_02);
fprintf('%-12.4f %-18.4f %-18.4f\n', 0.0100, Sd_01,   u_res_01);
fprintf('%-12.4f %-18.4f %-18.4f\n', 0.0020, Sd_002,  u_res_002);
fprintf('%-12.4f %-18.4f %-18.4f\n', 0.0005, Sd_0005, u_res_0005);

%% Plot
figure('Units','normalized','Position',[0.05 0.05 0.9 0.85]);

cases = {t_02,   u_02,   0.0200; ...
         t_01,   u_01,   0.0100; ...
         t_002,  u_002,  0.0020; ...
         t_0005, u_0005, 0.0005};

colors = {'#1f77b4','#d62728','#2ca02c','#9467bd'};

for j = 1:4
    t_j  = cases{j,1}';
    u_j  = cases{j,2};
    dt_j = cases{j,3};

    subplot(4,1,j);
    plot(t_j, u_j, 'Color', colors{j}, 'LineWidth', 1.2);
    ylabel('u (in)');
    title(sprintf('dt = %.4f s', dt_j));
    grid on; box on;
end

xlabel('Time (s)');
sgtitle('Problem 4: Displacement Time History', ...
    'FontSize', 13, 'FontWeight', 'bold');