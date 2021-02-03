InitialVelocity = zeros(3,1);
InitialPosition = zeros(3,1);


%% Model parameters:
% Vehicle parameters:
param.C_af = 133000;
param.C_ar = 98800;
param.m = 1650;
param.v = 22;
param.a = 1.11;
param.b = 1.59;
param.Iz = 2315.3;
% Sampling time:
param.dt = 0.1;
% Safety parameters:
% Lateral displacement y:
param.y_max = 0.9;
param.y_min = -0.9;
% Lateral velocity v:
param.v_max = 1.2;
param.v_min = -1.2;
% Yaw angle deviation Phi:
param.Phi_max = 0.05;
param.Phi_min = -0.05;
% Yaw rate dPhi:
param.dPhi_max = 0.3;
param.dPhi_min = -0.3;
% Input constraints (steering angle):
param.steer_max = 0.5; %pi/2;
param.steer_min = -0.5; %-pi/2;
% Disturbance range
param.rd_min = -0.02;
param.rd_max = 0.02;