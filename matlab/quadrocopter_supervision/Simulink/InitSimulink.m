%Initialization Script
RAD2DEG = 180 / pi;

% Sample times
%Tend = 35;
sim_dt = 0.0005;
%sim_dt = dt;

G_acc = 9.81;
Mass = 0.032;
InertiaMatrix = [
    16.571, 0.830, 0.718;
    0.831, 16.656, 1.800;
    0.718, 1.800, 29.262];
InertiaMatrix = InertiaMatrix * 1e-6;

Inv_InertiaMatrix = inv(InertiaMatrix);

body_arm = 0.04;

q0 = [1, 0, 0, 0]';
%q0 = eul2quat([pi/3, 0, 0], "XYZ");

Mix_XY = [
    -1, -1, 1,  1;
    -1,  1, 1, -1
    ];
Mix_Z = [1, -1, 1, -1];

%% Motors
motor_poly = [2.130295e-11 , 1.032633e-6 , 5.48456e-4];
motor_force_ub = 10;
motor_force2torque = 0.005964552;


%% Drags
linear_drag_coef = 5e-7;
quat_drag_coef = 5e-7;
rot_linear_drag_coef = 5e-7;

% linear_drag_coef = 0;
% quat_drag_coef = 0;
% rot_linear_drag_coef = 0;


%% Internal Controller
K_omega = 0.01;


%% Trajectory Object
%% 1) Computation of the path
% % Using trajectory class
% % Start from [-1.5, 1.0, 0] and go to [1.5, 0, 0] in 10 seconds.
% p0 = [0, 0.0, 0]';
% pe = [3.0, -1.0, 0]';
% Tend = 10;
% % Instantiate a Trajectory generator of degree 6.
% trajectory = TrajectoryClass(6);
% % Generate the polynomial that perform the interpolation. In this interface
% % I assume that I start from a still position and I end in a still
% % position (go to command).
% trajectory.generate(p0, pe, Tend);

coeff_x = trajectory.coeffs_px;
coeff_y = trajectory.coeffs_py;
coeff_z = trajectory.coeffs_pz;

M = der_mat(6)';
M2 = M^2;
M3 = M2 * M;
coeff_vx = M * coeff_x;
coeff_vy = M * coeff_y;
coeff_vz = M * coeff_z;

coeff_ax = M2 * coeff_x;
coeff_ay = M2 * coeff_y;
coeff_az = M2 * coeff_z;

coeff_jx = M3 * coeff_x;
coeff_jy = M3 * coeff_y;
coeff_jz = M3 * coeff_z;