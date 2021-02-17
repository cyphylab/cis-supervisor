% Import Experiment data
import_experimentdata;

time = experiment.time / 1e9;
N = length(time);

% Taking only position
CurrState = zeros(3, N);
ExpectedState = CurrState;
ActiveFlag = experiment.fieldactive;
AIndex = ActiveFlag > 0;
for (i = 1:N) 
    CurrState(1, i) = experiment.fieldstate_curr0(i);
    CurrState(2, i) = experiment.fieldstate_curr1(i);
    CurrState(3, i) = experiment.fieldstate_curr2(i);
    
    ExpectedState(1, i) = experiment.fieldstate_pred0(i);
    ExpectedState(2, i) = experiment.fieldstate_pred1(i);
    ExpectedState(3, i) = experiment.fieldstate_pred2(i);
end

%%
close all;
trj_plot = plot3(CurrState(1,:), CurrState(2,:), CurrState(3,:)); hold on;
act_plot = plot3(CurrState(1, AIndex), CurrState(2, AIndex), CurrState(3, AIndex), 'o');
legend('Trajectory', 'SupervisorActive');
set(trj_plot(1),'linewidth', 2);
set(act_plot(1),'linewidth', 2);
title('Trajectory');
xlabel('x [m]');
ylabel('y [m]');
zlabel('z [m]');

axis equal;
grid on;
