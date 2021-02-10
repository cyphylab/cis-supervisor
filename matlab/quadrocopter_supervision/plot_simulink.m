%% Plot Simulink output:
figure;
sim_output = out;

% Plot corrected trajectory:
plot3(sim_output.X_sim(:, 1), sim_output.X_sim(:, 2), sim_output.X_sim(:, 3),'LineWidth',3.5);
hold on;
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