%% Supervision and Simulation:
sim_output = SystemSimulation(DSets, mdl_dt, K_disc, x0, trajectory, 1 * Tend, dt, -1);

%% Plotting:
time_v = (0:sim_output.Nsim-1) * dt;

close all;

% Plot corrected trajectory:
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

% Plot obstacles and start/finish points:
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

% Plot nominal trajectory:
Nsim = Tend/dt;             % simulation steps
trj_pos = [];
for step = 1 : Nsim
    % Get trajectory:
    trj_pos = [trj_pos; (trajectory.trj_eval(step * dt, 0))'];
end
plot3(trj_pos(:, 1), trj_pos(:, 2), trj_pos(:, 3),'Color', [1 0.8 0], 'LineWidth',3.5);

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