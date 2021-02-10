function [sim_output] = SystemSimulation(DSets, mdl, K, x0, trajectory, Tend, dt, method)

Nsim = Tend/dt;             % simulation steps

Nsets = length(DSets);
x_curr = x0;

Ad = mdl.Ad;
Bd = mdl.Bd;
Nx = mdl.Nx;
Nu = mdl.Nu;

% Store trajectory details:
sim_output.X_sim = zeros(Nsim, Nx);
sim_output.E_sim = zeros(Nsim, Nx);
sim_output.U_nom = zeros(Nsim, Nu);
sim_output.U_corr = zeros(Nsim, Nu);
sim_output.opt_val = zeros(Nsim, 1);
sim_output.supervision_active = zeros(Nsim, 1);
sim_output.active_region = zeros(Nsim, 1);
sim_output.Nsim = Nsim;

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
    x_ref = [trj_pos; trj_vel; trj_acc];            % reference next state
    sim_output.X_sim(step, :) = x_curr;             % true current state
    err = x_ref - x_curr;                           % trajectory error
    sim_output.E_sim(step, :) = err;
    u_des = ctrl_fun(K, err) + trj_jrk;             % nominal input
    sim_output.U_nom(step, :) = u_des;
    
    % Check if nominal input is NaN
    if (isnan(u_des))
        disp(step);
        disp(Nsim);
        error("Nominal input is NaN!");
    end
    
    x_next_des = Ad * x_curr + Bd * u_des;    % nominal next state
    
    % Check if the nominal next state belongs to the CIS.
    [guard, ~] = isContained(x_next_des, DSets, []);
    % If desired next state belongs to the CIS we are good:
    if (guard)
        x_curr = x_next_des;
        sim_output.U_corr(step, :) = u_des;
        % Else we correct the nominal input:
    else
        sim_output.supervision_active(step) = step;
        
        % We solve for the union of CISs by checking if there exists a next
        % state in any of the CISs separately:
        u_cand = zeros(size(Bd,2), Nsets);  % candidate input
        f_cand = zeros(1, Nsets);           % candidate optimal cost
        for k = 1:Nsets
            % CIS inequalitites:
            %             cisA = DSets{k}.RCIS.A;
            %             cisb = DSets{k}.RCIS.b;
            %             Gu = DSets{k}.inputA;
            %             Fu = DSets{k}.inputb;
            %             % wrt to u:
            %             Aineq = [rcisA * Bd; Gu];
            %             bineq = [rcisb - rcisA * Ad * x_curr; Fu];
            
            % For the quadrocopter, the CIS inequalities wrt u
            % are simplified to box constraints:
            [lb,ub] = simplify2box(x_curr, DSets{k}, mdl);
            
            if (method==1)
                % Approach 1: use a solver:
                Aineq = [-eye(length(lb)); eye(length(ub))];
                bineq = [-lb; ub];
                % Cost: || u - udes ||^2.
                H = eye(size(Bd,2));
                c = -2 * u_des;
                % Call solver:
                result = solveGurobiQP(H, c, Aineq, bineq);
                
                if (strcmp(result.status,'OPTIMAL')==1)
                    u_cand(:, k) = result.x;
                    f_cand(k) = result.objval + norm(u_des)^2;
                else
                    u_cand(:, k) = NaN*ones(size(Aineq,2),1);
                    f_cand(k) = NaN;
                end
            else
                
                % Approach 2: Euclidean projection on a box has analytical
                % solution:
                [u_cand(:,k), f_cand(k)] = project2box(lb, ub, u_des);
            end
        end
        
        % Sanity check:
        if (min(isnan(f_cand)))
            disp(step);
            disp(Nsim);
            error("All the optimizations returned NaN..");
        end
        
        % Find minimum cost over all CISs:
        [sim_output.opt_val(step, :), next_cis] = min(f_cand);
        sim_output.active_region(step) = next_cis;
        % Select input corresponding to minimum cost:
        u_corrected = u_cand(:, next_cis);
        sim_output.U_corr(step, :) = u_corrected;
        
        % Make sure that next state is indeed in the next cis:
        x_next = Ad * x_curr + Bd * u_corrected;
        %         guard = DSets{next_cis}.Polyhedron.contains(x_next);
        if (~isContained(x_next, DSets, []))
            disp(step);
            disp(Nsim);
            cisA = DSets{next_cis}.CIS.A;
            cisb = DSets{next_cis}.CIS.b;
            sat = (cisA * x_next <= cisb);
            idcs = find(sat==0);
            diff = abs(cisA(idcs,:) * x_next - cisb(idcs));
            disp(diff)
            error("Next state not in a CIS!");
        end
        % Update state:
        x_curr = x_next;
    end
    
end


end