classdef IntegratorClass < handle
    properties
        state;
        size;
        input;
        AMat;
        BMat;
        dt;
    end
    
    methods
        function obj = IntegratorClass(Nx, dt)
            obj.state = zeros(Nx, 1);
            obj.input = 0;
            obj.size = Nx;
            obj.dt = dt;
            
            obj.genA(dt);
            obj.genB(dt);
        end
        
        function Init(obj, x0)
            obj.state(1) = x0;
        end
        
        function updateInput(obj, u)
            obj.input = u;
        end
        
        function genA(obj, dt)
            % Generate the Dynamics matrix A of the integrator system.
            intMatrix = eye(obj.size);
            fkt = 1;
            for (i = 1: obj.size - 1)
                fkt = fkt * i;
                intMatrix = intMatrix + diag(dt^i * ones(obj.size - i, 1), i)/fkt;
            end
            
            obj.AMat = intMatrix;
        end
        
        function genB(obj, dt)
            % Generate the Input matrix B of the integrator system.
            ctrlMatrix = zeros(obj.size, 1);
            ctrlMatrix(obj.size, 1) = dt;
            
            obj.BMat = ctrlMatrix;
        end
        
        function integrate(obj, dt)
            if (dt ~= obj.dt)
                obj.dt = dt;
                obj.genA(obj.dt);
                obj.genB(obj.dt);
            end
            obj.state = obj.AMat * obj.state +  obj.BMat * obj.input;
        end
        
        function A = getAMatrix(obj)
            A = obj.AMat;
        end
        
        function B = getBMatrix(obj)
            B = obj.BMat;
        end
        
    end
end
