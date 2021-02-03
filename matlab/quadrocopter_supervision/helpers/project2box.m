function [x_cand, f_cand] = project2box(lb, ub, x0)

x_cand = zeros(length(ub),1);

% Check feasibility:
if (min(ub-lb)<0)   % infeasible
    f_cand = NaN;
    x_cand = NaN*ones(length(ub),1);
else
    %           lb(j), if x0(j) < lb(j),
    % x*(k) =   ub(j), if x0(j) > ub(j),
    %           x0(j), otherwise.
    for j = 1:length(ub)
        if (x0(j)<lb(j))
            x_cand(j) = lb(j);
        elseif (x0(j)>ub(j))
            x_cand(j) = ub(j);
        else
            x_cand(j) = x0(j);
        end
    end
    f_cand = norm(x_cand-x0);
end

end