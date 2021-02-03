function [u_corr, opt_val, next_cis] = solveAnalytical(x0, udes, DSets, mdl)
Ad = mdl.Ad;
Bd = mdl.Bd;
Gu = DSets{1}.inputA;
Fu = DSets{1}.inputb;

u_cand = zeros(size(Gu,2),length(DSets));
f_cand = zeros(1,length(DSets));

for i = 1:length(DSets)
    
    cisA = DSets{i}.CIS.A;
    cisb = DSets{i}.CIS.b;
    
    Aineq = [cisA*Bd; Gu];
    bineq = [cisb - cisA*Ad*x0; Fu];
    
    % Normalize by the non-zero value of each row:
    for j = 1:size(Aineq,1)
        idx = find(Aineq(j,:)~=0);
        if (length(idx)>1)
            disp(i);
            disp(j);
            error("More than one non-zero value, not a hyper-rectangle");
        elseif (length(idx)==1)
            val = abs(Aineq(j,idx));
            Aineq(j,:) = Aineq(j,:)/val;
            bineq(j) = bineq(j)/val;
        else
            % Everything is zero, we skip.
        end
    end
    
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
    
    % Check feasibility:
    if (min(ub-lb)<0)   % infeasible
        f_cand(i) = Inf; 
    else
        for j = 1:size(Aineq,2)
            if (udes(j)<lb(j))
                u_cand(j,i) = lb(j);
            elseif (udes(j)>ub(j))
                u_cand(j,i) = ub(j);
            else
                u_cand(j,i) = udes(j);
            end
        end
        f_cand(i) = norm(u_cand(:,i)-udes);
    end
end

[opt_val, next_cis] = min(f_cand);
u_corr = u_cand(:,next_cis);

end