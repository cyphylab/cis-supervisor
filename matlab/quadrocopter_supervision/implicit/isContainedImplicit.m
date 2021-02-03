function [isContained, idx] = isContainedImplicit(x0, DSets, check_idx)

idx = zeros(length(DSets),1);
isContained = false;

if (isempty(check_idx))
    membership = zeros(length(DSets),1);
    parfor i = 1:length(DSets)
        Gx = DSets{i}.Gstate;
        Gu = DSets{i}.Ginput;
        Gv = DSets{i}.Gvirtual;
        F = DSets{i}.Fcl;
        % Feasibility problem:
        % minimize 0
        % s.t. Gu u + Gv v < F - Gx x0
        
        f = zeros(size(Gu,2)+size(Gv,2),1);
        A = [Gu Gv];
        b = F - Gx*x0;
        
        result = solveGurobiLP(f, A, b);
        
        if (strcmp(result.status,'OPTIMAL')==1)
            membership(i) = 1;
            idx(i) = 1;
        end
        
    end
    isContained = (max(membership)==1);
else
    Gx = DSets{check_idx}.Gstate;
    Gu = DSets{check_idx}.Ginput;
    Gv = DSets{check_idx}.Gvirtual;
    F = DSets{check_idx}.Fcl;
    % Feasibility problem:
    % minimize 0
    % s.t. Gu u + Gv v < F - Gx x0
    
    f = zeros(size(Gu,2)+size(Gv,2),1);
    A = [Gu Gv];
    b = F - Gx*x0;
    
    result = solveGurobiLP(f, A, b);
    
    if (strcmp(result.status,'OPTIMAL')==1)
        isContained = true;
        idx(check_idx) = 1;
    end
    
end