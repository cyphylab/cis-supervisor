Nsim = Tend/dt;             % simulation steps
x_curr = x0;
% Store trajectory details:
X_sim = zeros(Nsim, Nx);    % sequence of states
U_nom = zeros(Nsim, Nu);    % sequence of nominal inputs
U_corr = zeros(Nsim, Nu);   % sequence of corrected inputs
Opt_fval = zeros(Nsim, 1);  % optimization cost at each step
corr_idx = zeros(Nsim, 1);  % steps in which we correct
cis_idx = zeros(Nsim, 1);   % active CIS

% Let us solve the following:
% minimize {over u}     (u-u_des)^T (u-u_des)
% s.t.                  Gstate A x_curr + Gstate B u + Ginput u' + Gvirtual v < flifted
%                       Gu u < Fu
for step = 1 : Nsim
    % Get trajectory:
    trj_pos = trajectory.trj_eval(step * dt, 0);
	trj_vel = trajectory.trj_eval(step * dt, 1);
	trj_acc = trajectory.trj_eval(step * dt, 2);
    trj_jerk = trajectory.trj_eval(step * dt, 3);
    % Get nominal input:
    x_ref = [trj_pos; trj_vel; trj_acc];    % reference next state
    X_sim(step,:) = x_curr;                 % true current state
    err = x_ref - x_curr;                   % error
    u_des = ctrl_fun(K, err) + trj_jerk;    % nominal input
	U_nom(step, :) = u_des;
    % Check if nominal input is NaN
    if (isnan(u_des)) 
        disp(step);
        disp(Nsim);
        error("Nominal input NaN!");
    end
    x_next_des = Ad*x_curr + Bd*u_des;    % nominal next state
    % Check if the nominal next state belongs to the CIS:
    guard = checkContainmentImplicit(x_next_des, DSets,[]);
    
    % If desired next state belongs to the CIS we are good:
    if (guard == 1)
        x_curr = x_next_des;
    % Else we correct the nominal input:
    else
        corr_idx(step) = step;
        % We solve for the union of CISs by checking if there exists a next
        % state in any of the CISs separately:
        u_cand = zeros(size(Bd,2), Nsets);  % candidate input
        f_cand = zeros(1, Nsets);           % candidate optimal cost
        v_cand = zeros(size(Bd,2)*(L+1), Nsets);           % candidate virtual inputs
        for k = 1:Nsets            
            % Implicit CIS inequalitites:
            Gstate = DSets{k}.Gstate;
            Ginput = DSets{k}.Ginput;
            Gvirtual = DSets{k}.Gvirtual;
            Flifted = DSets{k}.flifted;
            % wrt to u:
%             Aineq = [Gstate*Bd Ginput Gvirtual; Gu zeros(size(Gu,1), size(Ginput,2)+size(Gvirtual,2))];
%             bineq = [Flifted - Gstate * Ad * x_curr; Fu];
            Aineq = [Gstate*Bd Ginput Gvirtual; blkdiag(Gu,Gu) zeros(12,size(Gvirtual,2))];
            bineq = [Flifted - Gstate * Ad * x_curr; Fu; Fu];
            
            % Original cost: || u - udes ||^2.
            HOG = blkdiag(eye(size(Gu,2)),zeros(size(Ginput,2)+size(Gvirtual,2)));
            cOG = [-2 * u_des; zeros(size(Ginput,2)+size(Gvirtual,2),1)];
            H = HOG;
            c = cOG;
            
            % Notice the formulation above, even though I have already
            % constraints on "the next" input (Ginput), I think that what I
            % want is to replace the current state with Ax+Bu and then
            % Ginput will refer to the input after u, and so on...

            % Call solver:
            result = solveGurobiQP(H, c, Aineq, bineq);

            if (strcmp(result.status,'OPTIMAL')==1)
                sol = result.x;
                fval = result.objval;
                u_cand(:, k) = sol(1:3);
                f_cand(k) = fval;
                v_cand(:,k) = sol(4:end);
            else
                u_cand(:, k) = [NaN; NaN; NaN];
                f_cand(k) = NaN;
                v_cand(:,k) = NaN*ones(3*L+3,1);
            end

        end
        
        % Sanity check:
        if (min(isnan(f_cand))) 
            disp(step);
            disp(Nsim);
            error("All the Optimizations returned NaN");
        end
        
        % Find minimum cost over all CISs:
        [Opt_fval(step, :), next_cis] = min(f_cand);
        cis_idx(step) = next_cis;
        % Select input corresponding to minimum cost:
        u_corrected = u_cand(:, next_cis);
        U_corr(step, :) = u_corrected;
        % Make sure that next state is indeed in the next cis:
        x_next = Ad * x_curr + Bd * u_corrected;
        guard = checkContainmentImplicit(x_next, DSets, next_cis);
        if (guard==0)
            disp(step);
            disp(Nsim);
            v_in = v_cand(:,next_cis);
            Gx = DSets{next_cis}.Gstate;
            Gi = DSets{next_cis}.Ginput;
            Gv = DSets{next_cis}.Gvirtual;
            Fl = DSets{next_cis}.flifted;
            sat = (Gx * x_next + [Gi Gv]*v_in <= Fl);
            idcs = find(sat==0);
            diff = abs(Gx(idcs,:) * x_next + [Gi(idcs,:) Gv(idcs,:)]*v_in - Fl(idcs));
            % check input constraints separately:
            min([blkdiag(Gu,Gu) zeros(12,size(Gv,2))]*[u_corrected; v_in] <= [Fu;Fu])
            error("Next state not in a CIS!");
        end
        % Update state:
        x_curr = x_next;
    end
    
end