function [Gu, Fu] = generateInputConstraints(mdl_cnstr)
Gu = [
    1  0  0; % x <  
   -1  0  0; % x >
    0  1  0; % y <
    0 -1  0; % y >
    0  0  1; % z <
    0  0 -1; % z >
    ];
Fu = mdl_cnstr.Jmax * ones(size(Gu,1),1);
end