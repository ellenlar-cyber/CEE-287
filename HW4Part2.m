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
tol = 0.0001; % 1% tolerance
% initialize Cy result vector
Cy = zeros(1,5);
%% Loop through ground motions
for i = 1:5
    ug = ground_motion{i};
    
    % Cy parameters
    Cy_min = 0;
    Cy_max = 1.0;
    Cy_step = -0.0001;

    % If no collapse, assume Cy(i) = Cy_max
    Cy(i) = Cy_max;
    
    % Create for loop to try Cy values
    for j = Cy_max: Cy_step: Cy_min
        Cy_try = j;

        [~, ~, ~, ~, Sd, ~] = Bilinear_SDOF_Response_NL(Tn, z, ug, dt, 0, 0, Cy_try, 'Cy', a, 'linear');

        m  = 1;
        wn = 2 * pi / Tn;
        k  = m * wn^2;
        c  = 2 * z * wn * m;
        g  = 981;             
    
        Fy = Cy_try * m * g;
        u_y = Fy / k;

        u_collapse = u_y * (1 - 1/a);

        if Sd >= u_collapse
            Cy1 = Cy_try; % Store the current Cy value
            break; % Exit the loop if the tolerance condition is met
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
%   initialize Ce result vector
    Ce = zeros(1,5);

% Loop through ground motions
for i = 1:5
    ug = ground_motion{i};
    % start with elastic condition as upper bound to iterate between
    Cy_elastic = 5;

    [~,~,~,~,Sd,~] = Bilinear_SDOF_Response_NL(Tn, z, ug, dt, 0, 0, Cy_elastic, 'Cy',a,'average');
    
    wn = 2 * pi / Tn;
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
% 
fprintf('Part (a): Minimum Cy\n')
fprintf('Cy_Parking = %.4f\n', Cy_parking)
fprintf('Cy_SLAC1   = %.4f\n', Cy_SLAC1)
fprintf('Cy_SLAC2   = %.4f\n', Cy_SLAC2)
fprintf('Cy_VA1     = %.4f\n', Cy_VA1)
fprintf('Cy_VA2     = %.4f\n', Cy_VA2)

fprintf('\nPart (b): Rc Values\n')
fprintf('Rc_Parking = %.4f\n', Rc_parking)
fprintf('Rc_SLAC1   = %.4f\n', Rc_SLAC1)
fprintf('Rc_SLAC2   = %.4f\n', Rc_SLAC2)
fprintf('Rc_VA1     = %.4f\n', Rc_VA1)
fprintf('Rc_VA2     = %.4f\n', Rc_VA2)
fprintf('Rc_avg     = %.4f\n', Rc_avg)    


    
