function [model_dt, model_constraints, model_ct] = InitModel3d(dt)
% Building the Model and Nominal Controller

%% Load the model parameters
model_constraints = vehicle_model();

%% Flat output modeling:
%%% (Since I always mess everything up with the indexes and orders...)
%%%
%%% Continuous-time system - 4th order integrator:
%%% This means that the higher derivative is the 4th and it's the input
%%% The state will be a vector of length 4 and will contain from the 0th to
%%% the 3rd derivative
%%%
%%% u -->|1/s|-- d3x/dt3 -->|1/s|-- d2x/dt2 -->|1/s|-- dx/dt -->|1/s|-- x
%%%
%%%       +-       -+
%%%       |   x     |
%%%       | dx/dt   |
%%%  X =  | d2x/dt2 |
%%%       | d3x/dt3 |
%%%       |         |
%%%       +-       -+
%%%
Nx = 3;

Ac_1D = diag(ones(Nx - 1, 1), 1);
Bc_1D = zeros(Nx, 1);
Bc_1D(Nx) = 1;
Cc_1D = eye(3);
Dc_1D = zeros(3,1);
fprintf('Creating ct integrator (Order: %d)\n', Nx);
disp('A matrix: ');
disp(Ac_1D);
disp('B Matrix: ');
disp(Bc_1D);
disp('');

ct_sys = ss(Ac_1D, Bc_1D, Cc_1D, Dc_1D);

% Discretization of the integrator
dt_sys = c2d(ct_sys, dt, 'zoh');
Ad_1D = dt_sys.A;
Bd_1D = dt_sys.B;
fprintf('Creating dt integrator (Order: %d, SampleTime: %f)\n', Nx, dt);
disp('A matrix: ');
disp(Ad_1D);
disp('B Matrix: ');
disp(Bd_1D);
disp('');



%% 3D Model
Nx_3D = 3 * Nx;
Nu = 3;

A_3D = zeros(Nx_3D);
A_3D(1:3, 4:6) = eye(3);
A_3D(4:6, 7:9) = eye(3);
B_3D = zeros(Nx_3D, Nu);
B_3D(7:9, :) = eye(3);


Ad_3D = eye(Nx_3D);
Ad_3D(1:3, 4:6) = eye(3) * Ad_1D(1,2);
Ad_3D(1:3, 7:9) = eye(3) * Ad_1D(1,3);

Ad_3D(4:6, 7:9) = eye(3) * Ad_1D(2,3);

Bd_3D = zeros(9, Nu);
Bd_3D(1:3, :) = eye(3) * Bd_1D(1);
Bd_3D(4:6, :) = eye(3) * Bd_1D(2);
Bd_3D(7:9, :) = eye(3) * Bd_1D(3);

model_ct = struct(...
    'A', A_3D, ...
    'B', B_3D, ...
    'A_1d', Ac_1D,...
    'B_1d', Bc_1D,...
    'Nx', Nx_3D, ...
    'Nu', Nu, 'dt', dt);


model_dt = struct(...
    'Ad', Ad_3D, ...
    'Bd', Bd_3D, ...
    'Ad_1d', Ad_1D,...
    'Bd_1d', Bd_1D,...
    'Nx', Nx_3D, ...
    'Nu', Nu, 'dt', dt);
end
