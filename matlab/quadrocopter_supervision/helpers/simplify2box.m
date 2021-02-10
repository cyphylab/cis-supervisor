function [lb, ub] = simplify2box(x0, DSet, mdl)
abs_tol = 1e-10; % absolute tolerance for checking containment.

Ad = mdl.Ad;
Bd = mdl.Bd;
rcisA = DSet.RCIS.A;
rcisb = DSet.RCIS.b;
Gu = DSet.inputA;
Fu = DSet.inputb;

Aineq = [rcisA*Bd; Gu];
bineq = [rcisb - rcisA*Ad*x0; Fu];

% Normalize by the non-zero value of each row:
guard = 0;
for j = 1:size(Aineq,1)
    idx = find(Aineq(j,:)~=0);
    if (length(idx)>1)
        disp(j);
        error("More than one non-zero value, not a hyper-rectangle");
    elseif (length(idx)==1)
        val = abs(Aineq(j,idx));
        Aineq(j,:) = Aineq(j,:)/val;
        bineq(j) = bineq(j)/val;
    else
        % Everything is zero, we check feasibility:
        if (bineq(j)<0 && abs(bineq(j))>abs_tol) % infeasible.
            guard = 1;
        end
    end
end

if (guard==1)
    ub = -ones(size(Aineq,2),1);
    lb = ones(size(Aineq,2),1);
else
    % Extract the box:
    lb = zeros(size(Aineq,2),1);
    ub = zeros(size(Aineq,2),1);
    for j = 1:size(Aineq,2)
        upper = (Aineq(:,j)>0);
        ubcol = bineq(upper);
        ub(j) = min(ubcol);
        
        lower = (Aineq(:,j)<0);
        lbcol = -bineq(lower);
        lb(j) = max(lbcol);
    end
end

end