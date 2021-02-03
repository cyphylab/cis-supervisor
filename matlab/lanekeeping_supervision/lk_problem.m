function [A,B,E,U,W,dyn,S] = lk_problem(param)
%% Define safe set for lateral-yaw model:

% All constraints are intervals:
Gstate = [eye(4); -eye(4)];
Fstate = [ param.y_max;   param.v_max;  param.Phi_max;  param.dPhi_max; 
          -param.y_min; -param.v_min; -param.Phi_min; -param.dPhi_min];
S = Polyhedron("H", [Gstate Fstate]);

%% Lateral-yaw dynamics (velocity is fixed). 
% x(t+1) = A x(t) + B u(t) + E d(t)

% Model parameters:
C_af = param.C_af;
C_ar = param.C_ar;
m = param.m;
v = param.v;
b = param.b;
Iz = param.Iz;
a = param.a;
dt = param.dt;

A = eye(4) + [  0   1                   v   0;
                0 -(C_af+C_ar)/m/v      0   (b*C_ar - a*C_af)/m/v-v;
                0   0                   0   1;
                0 (b*C_ar-a*C_af)/Iz/v  0  -(a^2*C_af + b^2*C_ar)/Iz/v]*dt;
            
B = [0; C_af/m; 0; a*C_af/Iz]*dt;

K = [0;0;0;0];

E = [0;0;-1;0]*dt;

% Input constraints (steering angle):
delta_f_max = param.steer_max;
delta_f_min = param.steer_min;
U = Polyhedron('H', [1 delta_f_max;
                     -1 -delta_f_min]);                 
XU = Polyhedron('H', [0  0  0  0  1  delta_f_max;
                      0  0  0  0 -1 -delta_f_min]);
                 
                 
% Disturbance constraints (desired yaw rate):
% d is the desired yaw rate, which is a time-varying external disturbance
% and computed from road curvature by d = vf /R0 where R0 is the signed
% radius of road curvature. Here we assume unmeasurable disturbance. 
rd_min = param.rd_min;
rd_max = param.rd_max;
W = Polyhedron('H', [1 rd_max; -1 -rd_min]);

dyn = Dyn(A, K, B, XU, [],[],[],{zeros(4,4)}, {E}, W);    