function result = solveGurobiLP(c, A, b)

model.obj = c;
model.sense = '<';
model.A = sparse(A);
model.rhs = b;
model.modelsense = 'min';
% Gurobi by default sets lb = 0, and ub = inf, so we have to set lb = -inf:
model.lb = -Inf*ones(size(A,2),1);
model.ub = Inf*ones(size(A,2),1);
% gurobi_write(model, 'gurobiSafeQP.lp');

params.FeasibilityTol = 1e-9;
% Default: 1e-6; Minimum: 1e-9; Maximum: 1e-2.
params.OptimalityTol = 1e-9;
% Dual feasibility tolerance: Default: 1e-6; Minimum: 1e-9; Maximum: 1e-2.
% Reduced costs must all be smaller than OptimalityTol in the improving 
% direction in order for a model to be declared optimal.
params.OutputFlag = 0;
result = gurobi(model,params);

% gurobi_feasrelax(model, relaxobjtype, minrelax, penalties, params)
%
% This function computes a feasibility relaxation for the input model 
% argument. The feasibility relaxation is a model that, when solved, 
% minimizes the amount by which the solution violates the bounds and 
% linear constraints of the original model. You must provide a penalty to 
% associate with relaxing each individual bound or constraint 
% (through the penalties argument). These penalties are interpreted in 
% different ways, depending on the value of the relaxobjtype argument.
% 
% 
% penalties: The penalties argument is a struct array:
%   lb Penalty for violating each lower bound.
%   ub Penalty for violating each upper bound.
%   rhs Penalty for violating each constraint.
% (default: all Inf).
% 
% relaxobjtype: The approach used to impose penalties on violations:
%   relaxobjtype=0, min the sum of the weighted magnitudes of violations.
%   relaxobjtype=1, min weighted sum of the squares of violations.
%   relaxobjtype=2, min the weighted count of violations.
%   Special penalty value Inf to indicate that the corresponding bound or 
%   constraint cannot be relaxed.
% 
% minrelax:
% 	minrelax=False, solution minimizes the cost of the violation. 
% 	minrelax=True, solution that minimizes the original objective, among 
%   solutions that minimize the cost of the violation (can be expensive).
% 
% EXAMPLE:
% 
%       penalties.rhs = ones(length(model.rhs),1);
%       feasrelaxresult = gurobi_feasrelax(model, 0, false, penalties);
% 
