
init;

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

%% Construct problem:
[Ad,Bd,Ed,U,W,dyn,S] = lk_problem(param);

safeA = S.A;
safeb = S.b;
Gu = U.A;
Fu = U.b;
Gw = W.A;
Fw = W.b;

%% Compute MCIS using PCIS
t = tic;
mcis = dyn.win_always(S,0,0,1);
tMCIS = toc(t);
volMCIS = mcis.volume;
disp(volMCIS);

% mprojYVdPhi = mcis.projection([1 2 4],'ifourier');
% mprojPhi = mcis.projection(3,'ifourier');

%% CIS2M:
L = 6;

% Project:
t = tic;
[cisA,cisb] = computeCIS(Ad, Bd, safeA, safeb, L, Gu, Fu, Ed, Gw, Fw);
cis = Polyhedron('H', [cisA cisb]);
tCIS = toc(t);
% disp(tCIS);
% cis = cis.minHRep;  % if computing all offline, remove redundant ineqs.
volCIS = cis.volume;
disp(volCIS);
% VolumePercentage = volCIS/volMCIS*100;
% disp(VolumePercentage)
% disp(cis.isEmptySet);

% projYVdPhi = cis.projection([1 2 4],'ifourier');
% projPhi = cis.projection(3,'ifourier');

%% Visualization
%
% figure(1);
% plot(mcis.slice(4,0),'alpha',0.5); hold on;
% % figure(2);
% plot(C.slice(4,0),'color','green')

% figure(1);
% plot(mcis.projection([1,2,3]),'alpha',0.5)
% % figure(2);
% hold on;
% plot(C.projection([1,2,3]),'color','green')
