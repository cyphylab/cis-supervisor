function [model] = vehicle_model()

% I consider a vehicle that can be controlled via angular rate
% and thrust derivative commands.

% Problem:
% The thrust derivative is annoying, because, in practice, has de bound 
% that depends on the thrust.
% It's a virtual control and I need to realized it using some trick. I was
% thinking that, being the discrete derivative df/dt = (f(k) - f(k-1)) / T,
% I can implement that input as:
% f(k) = df/dt * T + f(k-1)
% Indeed, in practice, I can control only the thrust.
%

% Pretending that this constraint does not exist, the constraint on the
% jerk are dependend only on the angular velocity.

% Mass
Mass = 0.032;
G = 9.81;

% Constraints
omega_max = 2000 / 180 * pi; % Gyro limits on each axis
T_max = 1.5 * G * Mass;      % Maximum collective thrust
T_min = 0.3 * G * Mass;      % Minimum collective thrust

% Maximum achievable speed and bound on each component
v_max = 1.5;  
v_max_i = v_max / sqrt(3);

% Considering that I can generate thrust in [0.3G, 1.5G] the acceleration
% along the vertical is in [-0.7G, 0.5G].
% This means that the sphere with the symmetric bounds is centered in zero
% with a radius of 0.5G.
a_max = (0.5 * G);
a_max_i = a_max / sqrt(3);


%% Computation of the constraints
jerk_max_i = (T_min / Mass) * omega_max / sqrt(3);

% Computation of the Jerk bound
% The norm of the jerk should not excede the (T_min * omega_max) and every
% component should not get more than 1/sqrt(3) * (T_min * omega_max).
% This is necessary to decouple the problem on separate axes without having
% the concern about the overall norm.
fprintf("Max jerk on each component: %3.3f\n", jerk_max_i);
fprintf("Max acceleration on each component: %3.3f\n", a_max_i);

model = struct('Mass', Mass,...
    'Vmax', v_max_i, 'Amax', a_max_i, 'Jmax', jerk_max_i,...
    'Tmax', T_max, 'Tmin', T_min);
end