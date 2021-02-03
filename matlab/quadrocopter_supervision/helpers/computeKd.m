function [Kd] = computeKd(K, dt)
    Kd = zeros(size(K));
    Kd(1) = K(1) + (2 * K(2) / dt) + (4 * K(3) / dt^2);
    Kd(2) = 2 * K(1) - 8 * K(3) / dt^2;
    Kd(3) = K(1) - (2 * K(2) / dt) + (4 * K(3) / dt^2);
end