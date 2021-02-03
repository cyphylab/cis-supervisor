%% Initialize:
init;           % clear variables, set paths. 
dt = 0.01;      % set sampling time.

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

% Create obstacle-free space:
obstacle_centers = [0.0 0.0 0.0; 1.5 1.0 0.0];
safe_distance = 0.5;
% Use IRIS:
% [S, Obs] = generateSafeRegion_IRIS(obstacle_centers, safe_distance, Gp_range, Fp_range);
% Semi-manual:
[Sp, Obs] = generateSafeRegion_Manual(obstacle_centers, safe_distance);

%% Compute the CIS:
% This should give a good CIS
% L = 500; dt = 0.001;

L = 15;
DSets = cell(length(Sp),1);
Ad = mdl_dt.Ad;
Bd = mdl_dt.Bd;
parfor i = 1:length(Sp)
    safeA = [G_state; [Sp(i).A zeros(size(Sp(i).A,1),6)]];
    safeb = [F_state; Sp(i).b];
    
    [Gstate,Ginput,Gvirtual,Fcl] = computeImplicitCIS(Ad, Bd, safeA, safeb, L, G_u, F_u);
    
    DSets{i} = struct('Index', i, ...
        'SafeSet', Polyhedron('H', [safeA, safeb]), ...
        'inputA', G_u, 'inputb', F_u, ...
        'Gstate', Gstate, 'Ginput', Ginput, 'Gvirtual', Gvirtual, 'Fcl', Fcl);
end
% % Check for CIS emptiness:
% for k = 1:length(DSets)
%     if (DSets{k}.CIS.isEmptySet)
%         fprintf('Empty set %d:', k);
%     end
% end
Nsets = length(DSets);
% DSetMap = generateDSetsMap(DSets);

% Check the projection of the velocity.

% Check if this point is contained:
%	I am interested in having some speed available
XX = zeros(9, 1);
XX(1:2:3) = -2.5;
XX(4) = 0.3;  % Speed 
isContainedImplicit(XX, DSets, [])
% 



%% Control Problem
%% 1) Computation of the path
% Using trajectory class
% Start from p_s with zero velocity and acceleration. 
% Arrive at p_f with zero velocity and acceleration. 
% In Tend seconds (??). 
p_s = [-2.0, -1.0, 0.0]';
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
guard = isContainedImplicit(x0, DSets, []);
if (guard==0)
    error('Point x0 not in the CIS');
end
guard = isContainedImplicit(xT, DSets, []);
if (guard==0)
    error('Point xT not in the CIS');
end

%% Controller design
% Piggy pole assignement...just have them negative.
gain = 1.5;
pCont = gain * [-2.0, -3.0, -4.0];
K_cont = place(mdl_ct.A_1d, mdl_ct.B_1d, pCont);

pDisc = exp(pCont.*dt); %poles in z-domain
K_disc = place(mdl_dt.Ad_1d, mdl_dt.Bd_1d, pDisc);

%% 2) Control
sim_output = SystemSimulationImplicit(L, DSets, mdl_dt, K_disc, x0, trajectory, 6 * Tend, dt);

system('say "Done ma frand!"')
%% Plotting:
time_v = (0:sim_output.Nsim-1) * dt;

close all;
figure;
index_correction = find(sim_output.supervision_active);
plot3(sim_output.X_sim(:, 1), sim_output.X_sim(:, 2), sim_output.X_sim(:, 3),'LineWidth',2.5); 
hold on;
plot3(sim_output.X_sim(index_correction, 1),...
    sim_output.X_sim(index_correction, 2),...
    sim_output.X_sim(index_correction, 3), 'ro', 'Color', 'cyan', 'LineWidth',2.5); 
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
hold off;

% figure;
% plot(sim_output.X_sim(:, 1), sim_output.X_sim(:, 2), 'LineWidth',2); hold on;
% cisses = {};
% counter = 1;
% for sim_step = 1 : 6
%     index_cis = (sim_output.supervision_active(sim_output.active_region == sim_step));
%     if (~isempty(index_cis))
%         col_ = [rand(), rand(), 0];
%         plot(...
%             sim_output.X_sim(index_cis, 1),...
%             sim_output.X_sim(index_cis, 2),...
%             'o', 'color', col_,'MarkerSize',10, 'LineWidth',1);
%         cisses(counter) = cellstr(['Supervision cis',mdl.Num2str(sim_step)]);
%         counter = counter + 1;
%     end
% end
% legend('Vehicle position', cisses{:}); 
% grid on;
% axis equal;
% xlabel("x");
% ylabel("y");
% zlabel("z");
% 
% figure;
% quiver(sim_output.X_sim(:, 1), sim_output.X_sim(:,2), U_nom(:,1), U_nom(:,2)); hold on;
% quiver(sim_output.X_sim(index_correction, 1), sim_output.X_sim(index_correction,2), U_corr(index_correction,1), U_corr(index_correction,2), 'r');
% grid on;
% axis equal;
% xlabel("x");
% ylabel("y");
% zlabel("z");
% legend("Vehicle position", "Supervision ON");
% 
figure;
for i = 1:3
    subplot(3, 1, i);
    plot(time_v, sim_output.U_nom(:, i),'LineWidth',1);
    title("Controls");
    hold on;
    plot(time_v, sim_output.U_corr(:, i),'-.','LineWidth',5);
    legend("Nominal", "Corrected");
    xlabel('time [s]');
    grid on;
end

for i = 1:3
    subplot(2, 3, i);
    plot(time_v, sim_output.X_sim(:, 3 + i),'LineWidth',1);
    title("Speed");
    xlabel('time [s]');
    grid on;
end

for i = 1:3
    subplot(2, 3, i + 3);
    plot(time_v, sim_output.X_sim(:, 6 + i),'LineWidth',1);
    title("Acceleration");
    xlabel('time [s]');
    grid on;
end

% 
% figure;
% plot(Opt_fval);
% legend("Optimization cost");

clear i