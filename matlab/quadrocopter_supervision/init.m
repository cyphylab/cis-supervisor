% Testing script.
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
addpath('../Library/cis2m_multi');
addpath('../Library/cis2m_multi/algorithms');
addpath('../Library/cis2m_multi/support_functions');