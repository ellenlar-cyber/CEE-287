%% HW 4 Number 2
clc; clear; close all
%% Load ground motions
data     = readmatrix('Ground_motions_assignment_1.csv');
t_record = data(:,1);
dt       = t_record(2) - t_record(1);
g        = 981;

% Columns: 1=time, 2=Parking, 3=SLAC-1, 4=SLAC-2, 5=VA-1, 6=VA-2
ag_Parking = data(:,2);
ag_SLAC1 = data(:,3);
ag_SLAC2 = data(:,4);
ag_VA1   = data(:,5);
ag_VA2   = data(:,6);

ground_motion = {ag_Parking,ag_SLAC1,ag_SLAC2,ag_VA1,ag_VA2};

%% Part a - Calculating minimum Cy
% Compute the minimum Seismic Coefficient bilinear SDOF system must have
    %given parameters:
   
    Tn = 1.0;           % seconds
    z = 0.05;           % damping ratio
    a = -0.15;

% Collapse will occur when Fs=0 (yield force is 0)
% Iterate through Cy values to find Fs = 0 (Cy min); 
    % similar iterative structure used in HW 2 part C

% Goal
Fs_goal = 0;
% Goal tolerance
Cy_tol = 0.0001; % 1% tolerance
% initialize Cy result vector
Cy = zeros(1,5);
%% Loop through ground motions
for i = 1:5
    ug = ground_motion{i};
    % Initialize variables for iterative test
    % start with elastic condition as upper bound to iterate between
    Cy1 = 5;
    % create an initial Cy at lower bound to iterate between
    Cy0 = 0.0000001;

    % while loop to continue iteration until tolerance acheived
    while abs(Cy1-Cy0) > Cy_tol
        Cy_1 = (Cy1+Cy0)/2;

        [~,ud,~, Fs,~,~] = Bilinear_SDOF_Response_NL(Tn, z, ug, dt, 0, 0, Cy_1, 'Cy',a,'average');
       
        % Change Cy1 based on result
        % if Fs is too high, Cy is too high, lower the upper bound
       % CORRECT - check if it DID collapse
    %% CHECK IF THIS IS THE PROPER ASSUMPTION
    % Fs hits zero AND displacement is still growing at that moment
       if any(abs(Fs(2:end)) < 0.001 & ud(2:end) > 0)
           Cy0 = Cy_1;  % collapsed → raise lower bound (need MORE strength)
       else
           Cy1 = Cy_1;  % survived → lower upper bound (can be WEAKER)
       end
    end
    Cy(i) = Cy1;
   
end

Cy_parking = Cy(1);
Cy_SLAC1 = Cy(2);
Cy_SLAC2 = Cy(3);
Cy_VA1 = Cy(4);
Cy_VA2 = Cy(5);

%% Part b - Corresponding maximum reduction factor Rc to avoid dynamic instability
% Rc is Ce (elastic)/ Cy(inelastic)
% Ce is Sa normalized by g

% Find Ce for all gm
  % initialize Ce result vector
    Ce = zeros(1,5);
%% Loop through ground motions
for i = 1:5
    ug = ground_motion{i};
    % start with elastic condition as upper bound to iterate between
    C = 5;

    [~,~,~,~,Sd,~] = Bilinear_SDOF_Response_NL(Tn, z, ug, dt, 0, 0, Cy_1, 'Cy',a,'average');
    wn = 2*pi/Tn;
    % Compute Sa from Sd
    Sa = wn^2 * Sd;
    % Normalize by g
    Ce(i) = Sa/g;
   
end

Ce_parking = Ce(1);
Ce_SLAC1 = Ce(2);
Ce_SLAC2 = Ce(3);
Ce_VA1 = Ce(4);
Ce_VA2 = Ce(5);

% Compute Rc Values
Rc_parking = Ce_parking/Cy_parking;
Rc_SLAC1 = Ce_SLAC1/Cy_SLAC1;
Rc_SLAC2 = Ce_SLAC2/Cy_SLAC2;
Rc_VA1 = Ce_VA1/Cy_VA1;
Rc_VA2 = Ce_VA2/Cy_VA2;

Rc_avg = mean([Rc_parking, Rc_SLAC1, Rc_SLAC2, Rc_VA1, Rc_VA2]);
    
fprintf('=== Part (a): Minimum Cy ===\n')
fprintf('Cy_Parking = %.4f\n', Cy_parking)
fprintf('Cy_SLAC1   = %.4f\n', Cy_SLAC1)
fprintf('Cy_SLAC2   = %.4f\n', Cy_SLAC2)
fprintf('Cy_VA1     = %.4f\n', Cy_VA1)
fprintf('Cy_VA2     = %.4f\n', Cy_VA2)

fprintf('\n=== Part (b): Rc Values ===\n')
fprintf('Rc_Parking = %.4f\n', Rc_parking)
fprintf('Rc_SLAC1   = %.4f\n', Rc_SLAC1)
fprintf('Rc_SLAC2   = %.4f\n', Rc_SLAC2)
fprintf('Rc_VA1     = %.4f\n', Rc_VA1)
fprintf('Rc_VA2     = %.4f\n', Rc_VA2)
fprintf('Rc_avg     = %.4f\n', Rc_avg)    
    
    
    
