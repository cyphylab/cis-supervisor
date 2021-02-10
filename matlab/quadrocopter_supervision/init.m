% Add required paths.
%
%
% Classes:
% - TrajectoryClass
% - ObstacleClass
% - Integrator3DClass
%
clear;
close all;
clc;

addpath('./helpers');
addpath('./classes');
addpath('./model');
addpath('./implicit');
addpath('./Simulink');
addpath('./Simulink/Functions');
addpath('./Simulink/Functions/Quaternion/');

addpath('../Library');
addpath('../Library/traj_helpers');
addpath('../../cis2m');
addpath('../../cis2m/algorithms');
addpath('../../cis2m/support_functions');