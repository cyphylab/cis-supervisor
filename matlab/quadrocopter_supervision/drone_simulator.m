%% Initialize:
init;           % clear variables, set paths. 
dt = 0.05; %125;      % set sampling time.

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

Ed = [1;1;1;0;0;0;0;0;0];
Gw = [1;-1];
Fw = [1e-3; 1e-3];
parfor i = 1:length(Sp)
    safeA = [G_state; [Sp(i).A zeros(size(Sp(i).A,1),6)]];
    safeb = [F_state; Sp(i).b];
    
%     [cisA, cisb] = computeCIS(mdl_dt.Ad, mdl_dt.Bd, safeA, safeb, T, L, G_u, F_u);
    [cisA, cisb] = computeCIS(mdl_dt.Ad, mdl_dt.Bd, safeA, safeb, T, L, G_u, F_u, Ed, Gw, Fw);
    cis = Polyhedron('H', [cisA, cisb]);
    cis = cis.minHRep;  % if computing all offline, remove redundant ineqs.
    
%     % Comment out for debugging:
%     UU = Polyhedron('H',[G_u F_u]);
%     if (~isInvariant(cis,UU,[],mdl_dt.Ad,mdl_dt.Bd))
%         warning('Result not numerically invariant.');
%     else
%         disp('Invariant in original space!')
%     end
    
    DSets{i} = struct('Index', i, ...
        'SafeSet', Polyhedron('H', [safeA, safeb]), ...
        'inputA', G_u, 'inputb', F_u, ...
        'CIS', cis);    
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
% In Tend seconds (??). 
p_s = [-1.0, -1.0, 0.0]';
p_f = [2.5, 2.0, 0.45]';
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

%% Plot Simulink output:
figure;
sim_output = out;
plot3(sim_output.X_sim(:, 1), sim_output.X_sim(:, 2), sim_output.X_sim(:, 3),'LineWidth',3.5);
hold on;
grid on;
axis equal;
xlabel("x");
ylabel("y");
zlabel("z");
for i = 1:length(Obs)
    plot(Obs(i),'color',[0.6350 0.0780 0.1840],'alpha',0.3);
    cisA = Obs(i).chebyCenter();
    plot3(cisA.x(1),cisA.x(2),cisA.x(3),'d','Color','red','LineWidth',5);
end
plot3(x0(1),x0(2),x0(3),'s','Color','magenta','LineWidth',2.5);
text(x0(1),x0(2),x0(3),'\leftarrow start','FontSize',14);
plot3(xT(1),xT(2),xT(3),'o','Color','green','LineWidth',2.5);
text(xT(1),xT(2),xT(3),'\leftarrow finish','FontSize',14);
xlim([-3,3]);
ylim([-3,3]);

Nsim = Tend/dt;             % simulation steps
trj_pos = [];
for step = 1 : Nsim
    % Get trajectory:
    trj_pos = [trj_pos; (trajectory.trj_eval(step * dt, 0))'];
end
plot3(trj_pos(:, 1), trj_pos(:, 2), trj_pos(:, 3),'Color', [1 0.8 0], 'LineWidth',3.5);

hold off;
