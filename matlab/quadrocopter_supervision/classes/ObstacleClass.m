classdef ObstacleClass < handle
    properties
        Tr = struct('R', [], 't', []);
    end
    
    methods
        function obj = ObstacleClass()
            obj.Tr.R = eye(3);
            obj.Tr.t = zeros(3,1);
        end
        
        function obj = setTransform(obj, Rotation, Translation)
            % Set the position and orientation of the Reference Frame
            % attached to the obstacle.
            obj.Tr.R = Rotation;
            obj.Tr.t = reshape(Translation, 3, 1);
        end
        
        function out = trj_eval(obj, val)
            % Evaluate the 'der' derivative of the trajectory at time 't'
            % in the reference 'ref'.
            % The 'ref' can be 'World' or 'Obstacle'.
            % Note: The 0 derivative is just the trajectory.
            if (der == 0)
                out = obj.Tr.R * (val - obj.Tr.t);
            else
                out = obj.Tr.R * val;
            end
            
        end
        
        function [x_rad, u_rad, T] = radial_eval(obj, x, u)
            % Compute the "Radial Space" form of the 3D state 'x' and 3D
            % input 'u'.
            %
            % x = Vector of size (3, StateLength);
            % u  = Vector of size 3
            Nx = size(x, 2);
            
            % Preallocate output vector
            x_rad = zeros(Nx, 1);
            u_rad = 0;  
            T = 0;
            
            % 3D Position of the vehicle with respect to the obstacle
            p_0 = obj.Tr.R * (x(:, 1) - obj.Tr.t);
            
            for ind = 1 : Nx
                % Radial Distance
                if (ind == 1)
                    x_rad(ind) = norm(p_0);
                end
                
                % Radial Velocity
                if (ind == 2)
                    p_1 = obj.Tr.R * x(:, 2); % First derivative
                    x_rad(ind) = dot(p_0, p_1)/x_rad(1);
                end
                
                % Radial Acceleration
                if (ind == 3)
                    p_1 = obj.Tr.R * x(:, 2); % 1st Derivative
                    p_2 = obj.Tr.R * x(:, 3); % 2nd Derivative
                    x_rad(ind) = (norm(p_1)^2 + dot(p_0, p_2))/x_rad(1) - ...
                        dot(p_0, p_1)/(x_rad(1)^2)*x_rad(2);
                end
                
                % Radial Jerk
                if (ind == 4)
                    p_1 = obj.Tr.R * x(:, 2); % 1st Derivative
                    p_2 = obj.Tr.R * x(:, 3); % 2nd Derivative
                    p_3 = obj.Tr.R * x(:, 4); % 3nd Derivative
                    x_rad(ind) = 3 * dot(p_1, p_2)/x_rad(1) + ...
                        dot(p_0, p_3) / x_rad(1) - ...
                        2.0 * (norm(p_1)^2 + dot(p_0, p_2)) * x_rad(2) / x_rad(1)^2 + ...
                        2.0 * dot(p_0, p_1) / x_rad(1)^3 * x_rad(2)^2 - ...
                        dot(p_0, p_1) / x_rad(1)^2 * x_rad(3);
                end
                
                % Radial Snap
                if (ind == 5)
                    disp('Something is rotten in LA')
                    p_1 = obj.Tr.R * x(:, 2); % 1st Derivative
                    p_2 = obj.Tr.R * x(:, 3); % 2nd Derivative
                    p_3 = obj.Tr.R * x(:, 4); % 3nd Derivative
                    p_4 = obj.Tr.R * x(:, 5); % 4th Derivative
                    x_rad(ind) = (4 * dot(p_1, p_3) + 3 * (dot(p_2, p_2)) + dot(p_0, p_4))/x_rad(1) - ...
                        3 * (3 * dot(p_1, p_2) + dot(p_0, p_3))/x_rad(1)^2 * x_rad(2) - ...
                        3 * (dot(p_1, p_1) + dot(p_0, p_2))/x_rad(1)^2 * x_rad(3) - dot(p_0, p_1)/x_rad(1)^2 * x_rad(4) + ...
                        6 * (dot(p_1, p_1) + dot(p_0, p_2))/x_rad(1)^3 * x_rad(2)^2 + 6 * dot(p_0, p_1)/x_rad(1)^3 * x_rad(3) * x_rad(2) - ...
                        6 * dot(p_0, p_1)/x_rad(1)^4 * x_rad(2)^3;
                end
            end
            
            % Input in the obstacle frame
            u_obst = obj.Tr.R * u;
            u_part = dot(p_0, u_obst) / x_rad(1);
            
            switch (Nx)
                case 1
                    T = 0;
                    u_rad = u_part + T;
                case 2
                    T = dot(p_1, p_1) / x_rad(1) - ...
                        dot(p_0, p_1) / (x_rad(1)^2) * x_rad(2);
                    u_rad = u_part + T;
                case 3
                    T = 3 * dot(p_1, p_2)/x_rad(1) - ...
                        2.0 * (norm(p_1)^2 + dot(p_0, p_2)) * x_rad(2) / x_rad(1)^2 + ...
                        2.0 * dot(p_0, p_1) / x_rad(1)^3 * x_rad(2)^2 - ...
                        dot(p_0, p_1) / x_rad(1)^2 * x_rad(3);
                    u_rad = u_part + T;
                case 4
                    T = (4 * dot(p_1, p_3) + 3 * (dot(p_2, p_2)))/x_rad(1) - ...
                        3 * (3 * dot(p_1, p_2) + dot(p_0, p_3))/x_rad(1)^2 * x_rad(2) - ...
                        3 * (dot(p_1, p_1) + dot(p_0, p_2))/x_rad(1)^2 * x_rad(3) - dot(p_0, p_1)/x_rad(1)^2 * x_rad(4) + ...
                        6 * (dot(p_1, p_1) + dot(p_0, p_2))/x_rad(1)^3 * x_rad(2)^2 + 6 * dot(p_0, p_1)/x_rad(1)^3 * x_rad(3) * x_rad(2) - ...
                        6 * dot(p_0, p_1)/x_rad(1)^4 * x_rad(2)^3;
                    u_rad = u_part + T;
                otherwise 
                    disp('Something is rotten in LA')
            end
        end
        
        
        function r = get_r(obj, pos)
            r = obj.Tr.R * (pos - obj.Tr.t);
        end
        
    end
end
