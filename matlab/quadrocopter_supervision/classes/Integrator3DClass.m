classdef Integrator3DClass < handle
    properties
		state;
		size;
		input;
		Ints;
    end
    
    methods
        function obj = Integrator3DClass(Nx, dt)
            % Constructor: Creates a Integrator3dClass with a state of
            % length Nx on each dimension.
            obj.state = zeros(3, Nx);
			obj.input = zeros(3, 1);
	    	obj.size = Nx;
			obj.Ints = [IntegratorClass(Nx, dt), IntegratorClass(Nx, dt), IntegratorClass(Nx, dt)]';
        end
        
        function Init(obj, xyz0)
            % Initialize the position of the integrator
			for (i = 1: 3)
				obj.Ints(i).Init(xyz0(i));
				obj.state(i, :) = obj.Ints(i).state;
			end
        end

		function updateInput(obj, u)
            % Update the input signal of the integrator
			obj.input = reshape(u, 3, 1);
			for (i = 1:3)
				obj.Ints(i).updateInput(obj.input(i))
			end
		end

		function integrate(obj, dt)
            % Integration step with a time step = dt
			for (i = 1:3)
				obj.Ints(i).integrate(dt);
				obj.state(i,:) = obj.Ints(i).state;
			end
        end

    end
end
