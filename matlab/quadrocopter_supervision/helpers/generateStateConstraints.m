function [G_state, F_state] = generateStateConstraints(mdl_cnstr, space_range)
% Velocity constraints:
Gv = [
    1  0  0; % x <  
   -1  0  0; % x >
    0  1  0; % y <
    0 -1  0; % y >
    0  0  1; % z <
    0  0 -1; % z >
    ];
Fv = mdl_cnstr.Vmax * ones(size(Gv,1),1);
% Acceleration constraints:
Ga = [
    1  0  0; % x <  
   -1  0  0; % x >
    0  1  0; % y <
    0 -1  0; % y >
    0  0  1; % z <
    0  0 -1; % z >
    ];
Fa = mdl_cnstr.Amax * ones(size(Ga,1),1);
% Operational volume in position space:

Gp_range = [
    1  0  0; % x <
   -1  0  0; % x >
    0  1  0; % y <
    0 -1  0; % y >
    0  0  1; % z <
    0  0 -1; % z >
    ];
Fp_range = space_range * ones(size(Gp_range,1),1);
% Concatenate state constraints:
G_state = blkdiag(Gp_range, Gv, Ga);
F_state = [Fp_range; Fv; Fa];
end
