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
% s.t.                  cisA ( A x_curr + B u ) < cisb
%                       Gu u < Fu
for step = 1 : Nsim
    % Get trajectory:
    trj_pos = trajectory.trj_eval(step * dt, 0);
    trj_vel = trajectory.trj_eval(step * dt, 1);
    trj_acc = trajectory.trj_eval(step * dt, 2);
    trj_jrk = trajectory.trj_eval(step * dt, 3);
    % Get nominal input:
    x_ref = [trj_pos; trj_vel; trj_acc];    % reference next state
    X_sim(step,:) = x_curr;                 % true current state
    err = x_ref - x_curr;                   % error
    u_des = ctrl_fun(K, err) + trj_jrk;     % nominal input
    U_nom(step, :) = u_des;
    % Check if nominal input is NaN
    if (isnan(u_des))
        disp(step);
        disp(Nsim);
        error("Nominal input NaN!");
    end
    x_next_des = Ad*x_curr + Bd*u_des;    % nominal next state
    % Check if the nominal next state belongs to the CIS:
    guard = 0;
    for k = 1:length(DSets)
        if (DSets{k}.CIS.contains(x_next_des))
            guard = 1;
        end
    end
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
        for k = 1:Nsets
            if (DSets{k}.CIS.contains(x_curr))
                % CIS inequalitites:
                cisA = DSets{k}.CIS.A;
                cisb = DSets{k}.CIS.b;
                % wrt to u:
                Aineq = [cisA * Bd; Gu];
                bineq = [cisb - cisA * Ad * x_curr; Fu];
                % Original cost: || u - udes ||^2.
                HOG = eye(size(Bd,2));
                cOG = -2 * u_des;
                H = HOG;
                c = cOG;
                
                % Call solver:
                result = solveGurobiQP(H, c, Aineq, bineq);
                
                if (strcmp(result.status,'OPTIMAL')==1)
                    sol = result.x;
                    fval = result.objval;
                    u_cand(:, k) = sol(1:3);
                    f_cand(k) = fval;
                else
                    u_cand(:, k) = [NaN; NaN; NaN];
                    f_cand(k) = NaN;
                end
            else
                u_cand(:, k) = [NaN; NaN; NaN];
                f_cand(k) = NaN;
            end
        end
        
        % Sanity check:
        if (min(isnan(f_cand)))
            disp(step);
            disp(Nsim);
            error("All the optimizations returned NaN..");
        end
        
        % Find minimum cost over all CISs:
        [Opt_fval(step, :), next_cis] = min(f_cand);
        cis_idx(step) = next_cis;
        % Select input corresponding to minimum cost:
        u_corrected = u_cand(:, next_cis);
        U_corr(step, :) = u_corrected;
        % Make sure that next state is indeed in the next cis:
        x_next = Ad * x_curr + Bd * u_corrected;
        %         guard = DSets{next_cis}.Polyhedron.contains(x_next);
        guard = 0;
        for k = 1:length(DSets)
            if (DSets{k}.CIS.contains(Ad * x_curr + Bd * u_corrected))
                guard = 1;
            end
        end
        if (guard==0)
            disp(step);
            disp(Nsim);
            %             cisA = DSets{next_cis}.CIS.A;
            %             cisb = DSets{next_cis}.CIS.b;
            %             sat = (cisA * x_next <= cisb);
            %             idcs = find(sat==0);
            %             diff = abs(cisA(idcs,:) * x_next - cisb(idcs));
            error("Next state not in a CIS!");
        end
        % Update state:
        x_curr = x_next;
        
        %         Nsim = Tend/dt;
        trajectory.generate(x_curr, xT, Tend); % - step*dt);
    end
    
end