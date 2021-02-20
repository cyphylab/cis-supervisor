% Import Experiment data
import_experimentdata;
import_gyro;


%%
%{
The data from the experiment contains:
Commanded jerk in body frame body: jerk_body
Commanded thrust: thrust
Commanded yaw rate: yaw_rate
Commanded angular velocity in body frame: ang_velocity
Current flat output state of the drone: state_curr 
Predicted flat output state of the drone: state_pred 
Active Supervision Flag: active

Matlab imports those filed as "field<nameofthefield>".
%}

time = experiment.time / 1e9;
time_meas = gyromeas.time / 1e9;

t00 = min(time(1), time_meas(1));
time = time - t00;
time_meas = time_meas - t00;

N = length(time);
Nmeas = length(time_meas);

CurrState = zeros(9, N);
GyroMeas = zeros(3, Nmeas);

CmdGyro = zeros(3, N);

ExpectedState = CurrState;
ActiveFlag = experiment.fieldactive;

UJerk = zeros(3, N);


AIndex = ActiveFlag > 0;
for (i = 1:N) 
    CurrState(1, i) = experiment.fieldstate_curr0(i);
    CurrState(2, i) = experiment.fieldstate_curr1(i);
    CurrState(3, i) = experiment.fieldstate_curr2(i);
    CurrState(4, i) = experiment.fieldstate_curr3(i);
    CurrState(5, i) = experiment.fieldstate_curr4(i);
    CurrState(6, i) = experiment.fieldstate_curr5(i);
    CurrState(7, i) = experiment.fieldstate_curr6(i);
    CurrState(8, i) = experiment.fieldstate_curr7(i);
    CurrState(9, i) = experiment.fieldstate_curr8(i);
    
    
    ExpectedState(1, i) = experiment.fieldstate_pred0(i);
    ExpectedState(2, i) = experiment.fieldstate_pred1(i);
    ExpectedState(3, i) = experiment.fieldstate_pred2(i);
    ExpectedState(4, i) = experiment.fieldstate_pred3(i);
    ExpectedState(5, i) = experiment.fieldstate_pred4(i);
    ExpectedState(6, i) = experiment.fieldstate_pred5(i);
    ExpectedState(7, i) = experiment.fieldstate_pred6(i);
    ExpectedState(8, i) = experiment.fieldstate_pred7(i);
    ExpectedState(9, i) = experiment.fieldstate_pred8(i);
    
    UJerk(1, i) = experiment.fieldjerk_body0(i);
    UJerk(2, i) = experiment.fieldjerk_body1(i);
    UJerk(3, i) = experiment.fieldjerk_body2(i);
    
    
end

CmdGyro(1, :) = experiment.fieldang_velocity0;
CmdGyro(2, :) = -experiment.fieldang_velocity1;
CmdGyro(3, :) = experiment.fieldang_velocity2;
    
GyroMeas(1, :) = gyromeas.fieldvalues0;
GyroMeas(2, :) = gyromeas.fieldvalues1;
GyroMeas(3, :) = gyromeas.fieldvalues2;


%% Check on the estimator
CheckEstVel = CurrState(4:6, :);
for (i = 1:N-1)
    dt = time(i+1) - time(i);
    CheckEstVel(:, i+1) = CurrState(4:6, i) + CurrState(7:9, i) * dt;
    % XXX I should save the quaternion or the Jerk in Inertial Frame....
end


%%
close all;
trj_plot = plot3(CurrState(1,:), CurrState(2,:), CurrState(3,:)); hold on;
act_plot = plot3(CurrState(1, AIndex), CurrState(2, AIndex), CurrState(3, AIndex), 'o');

% Plot obstacles and start/finish points:
%{
for i = 1:length(Obs)
    plot(Obs(i),'color',[0.6350 0.0780 0.1840],'alpha',0.3);
    cisA = Obs(i).chebyCenter();
    plot3(cisA.x(1),cisA.x(2),cisA.x(3),'d','Color','red','LineWidth',5);
end
%}
legend('Trajectory', 'SupervisorActive');
set(trj_plot(1),'linewidth', 2);
set(act_plot(1),'linewidth', 2);
title('Trajectory');
xlabel('x [m]');
ylabel('y [m]');
zlabel('z [m]');
axis equal;
grid on;

%%

figure;
for (i = 1 : 3)
    subplot(3,1,i);
    velocity_plot = plot(CurrState(i, :), 'o');
    hold on; plot(ExpectedState(i, :));
    plot(find(AIndex), CurrState(i, AIndex), 'og', 'Markersize', 12);
    title('Position Estimate');
    grid on;
    legend('curr', 'pred', 'supervision');
end

figure;
for (i = 1 : 3)
    subplot(1,3,i);
    velocity_plot = plot(CurrState(3 + i, :), 'o');
    hold on; plot(ExpectedState(3 + i, :));
    plot(find(AIndex), CurrState(3 + i, AIndex), 'og', 'Markersize', 12);
    title('Velocity Estimate');
    grid on;
    legend('curr', 'pred', 'supervision');
end

figure;
for (i = 1 : 3)
    subplot(1,3,i);
    velocity_plot = plot(CurrState(6 + i, :), 'o');
    hold on; plot(ExpectedState(6 + i, :));
    plot(find(AIndex), CurrState(6 + i, AIndex), 'og', 'Markersize', 12);
    title('Acceleration Estimate');
    grid on;
    legend('curr', 'pred', 'supervision');
end

figure;
for (i = 1 : 3)
    subplot(3,1,i);
    jerk_plot = plot(UJerk(i, :)); hold on;
    plot(find(AIndex), UJerk(i, AIndex), 'og', 'Markersize', 12);
    title('Jerk Command');
    grid on;
    legend('command', 'supervision');
end

figure;
for (i = 1 : 3)
    subplot(3,1,i);
    jerk_plot = plot(CurrState(3 + i, :)); hold on;
    plot(CheckEstVel(i, :)); hold on;
    plot(find(AIndex), CurrState(3 + i, AIndex), 'og', 'Markersize', 12);
 
    title('Check the estimation');
    grid on;
    legend('Estimated', 'A posteriori', 'supervision');
end

figure;
for (i = 1 : 2)
    subplot(2, 1, i);
    gyro_plot = stairs(time, 180 * CmdGyro(i, :) / pi); hold on;
    plot(time_meas, GyroMeas(i, :));
    plot(time(find(AIndex)), CmdGyro(i, AIndex), 'og', 'Markersize', 12);
    title('Gyros');
    grid on;
    legend('Commanded', 'Measured', 'supervision');
end
