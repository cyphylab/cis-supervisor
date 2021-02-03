function [isContained, idx] = checkContainmentImplicit(x0, DSets, check_idx)

idx = [];
isContained = false;

if (isempty(check_idx))
    for i = 1:length(DSets)
        Gx = DSets{i}.Gstate;
        Gu = DSets{i}.Ginput;
        Gv = DSets{i}.Gvirtual;
        F = DSets{i}.flifted;
        % Feasibility problem:
        % minimize 0
        % s.t. Gu u + Gv v < F - Gx x0
        
        f = zeros(size(Gu,2)+size(Gv,2),1);
        A = [Gu Gv];
        b = F - Gx*x0;
        
        result = solveGurobiLP(f, A, b);
        
        if (strcmp(result.status,'OPTIMAL')==1)
            isContained = true;
            idx = [idx; i];
        end
        
    end
else
    Gx = DSets{check_idx}.Gstate;
    Gu = DSets{check_idx}.Ginput;
    Gv = DSets{check_idx}.Gvirtual;
    F = DSets{check_idx}.flifted;
    % Feasibility problem:
    % minimize 0
    % s.t. Gu u + Gv v < F - Gx x0
    
    f = zeros(size(Gu,2)+size(Gv,2),1);
    A = [Gu Gv];
    b = F - Gx*x0;
    
    result = solveGurobiLP(f, A, b);
    
    if (strcmp(result.status,'OPTIMAL')==1)
        isContained = true;
    end
    
end