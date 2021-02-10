%% Initialize:
init;               % clear variables, set paths. 
dt = 0.10; %125;    % set sampling time.

% Define vehicle constraints and return the flat output model.
% mdl_dt: discrete time model
% mdl_cnstr: constraints of the vehicle
% mdl_ct: continuous time model (in case I need)
[mdl_dt, mdl_cnstr, mdl_ct] = InitModel3d(dt);

%% Define safe set:
space_range = 3.0;
% State Constraints:
[G_state, F_state] = generateStateConstraints(mdl_cnstr, space_range);
% Input constraints:
[G_u, F_u] = generateInputConstraints(mdl_cnstr);

% Create obstacle-free space:d
obstacle_centers = [0.0 0.0 0.0; 1.5 1.0 0.0];
safe_distance = 0.5;
% Use IRIS:
% [S, Obs] = generateSafeRegion_IRIS(obstacle_centers, safe_distance, Gp_range, Fp_range);
% Semi-manual:
[Sp, Obs] = generateSafeRegion_Manual(obstacle_centers, safe_distance);

%% Compute the CISs:
L = 6; T = 0;
DSets = cell(length(Sp),1);

% Define disturbance:
k = mdl_dt.Nx;
Ed = eye(k);
Gw = [eye(k); -eye(k)];

Fw_pos = 1e-2 * ones(3,1);
Fw_vel = 1e-2 * ones(3,1);
Fw_acc = 5e-2 * ones(3,1);
Fw_pva = [Fw_pos; Fw_vel; Fw_acc];

Fw = [Fw_pva; Fw_pva];
W = Polyhedron('H', [Gw Fw]);
parfor i = 1:length(Sp)
    safeA = [G_state; [Sp(i).A zeros(size(Sp(i).A,1),6)]];
    safeb = [F_state; Sp(i).b];
    
%     [cisA, cisb] = computeCIS(mdl_dt.Ad, mdl_dt.Bd, safeA, safeb, T, L, G_u, F_u);
    [cisA, cisb] = computeCIS(mdl_dt.Ad, mdl_dt.Bd, safeA, safeb, T, L, G_u, F_u, Ed, Gw, Fw);
    cis = Polyhedron('H', [cisA, cisb]);
    cis = cis.minHRep;  % if computing all offline, remove redundant ineqs.
    
    % This is to be used for the optimization problem:
    % Ax + Bu + E W \subseteq CIS <-> Ax + Bu \subseteq CIS - E W
    % <-> Ax + Bu \in RCIS.
    rcis = cis - Ed * W;
    
    DSets{i} = struct('Index', i, ...
        'SafeSet', Polyhedron('H', [safeA, safeb]), ...
        'inputA', G_u, 'inputb', F_u, ...
        'CIS', cis, 'RCIS', rcis);    
end
% Check for CIS emptiness:
for k = 1:length(DSets)
    if (DSets{k}.CIS.isEmptySet)
        fprintf('Empty set: %d. \n', k);
    end
end
% Clear temporary variables
clear cis cisA cisb safeA safeb i k

%% Computation of the path
% Using trajectory class
% Start from p_s with zero velocity and acceleration. 
% Arrive at p_f with zero velocity and acceleration. 
% In Tend seconds. 
p_s = [-1.5, -1.0, 0.0]';
p_f = [2.5, 2.0, 0.0]';
Tend = 20;
% Instantiate a Trajectory generator of degree 6.
trajectory = TrajectoryClass(6);
% Generate the polynomial that performs the interpolation. 
trajectory.generate(p_s, p_f, Tend);
% Start state (p_s,0,0):
x0 = zeros(mdl_dt.Nx, 1);
x0(1:3) = p_s;
% End state (p_f,0,0):
xT = zeros(mdl_dt.Nx, 1);
xT(1:3) = p_f;

% Check if x0, xT are in the CIS:
if (~isContained(x0, DSets, []))
    error('Point x0 not in the CIS');
end

if (~isContained(xT, DSets, []))
    error('Point xT not in the CIS');
end

%% Controller design:
% Piggy pole assignement...just have them negative.
gain = 1.5;
ctrl_dt = dt;
pCont = gain * [-1.0, -1.3, -1.2];
K_cont = place(mdl_ct.A_1d, mdl_ct.B_1d, pCont);

pDisc = exp(pCont.*ctrl_dt); %poles in z-domain
K_disc = place(mdl_dt.Ad_1d, mdl_dt.Bd_1d, pDisc);


Plant_dt = ss(mdl_dt.Ad_1d, mdl_dt.Bd_1d, eye(3), []);

%% Initialize the Simulink environment:
InitSimulink;

export_data;

