function [sim_output] = SystemSimulationImplicit(L, DSets, mdl, K, x0, trajectory, Tend, dt)

Nsim = Tend/dt;             % simulation steps
Nsets = length(DSets);
x_curr = x0;
Ad = mdl.Ad;
Bd = mdl.Bd;
n = mdl.Nx;
m = mdl.Nu;

% Store trajectory details:
sim_output.X_sim = zeros(Nsim, n);
sim_output.E_sim = zeros(Nsim, n);
sim_output.U_nom = zeros(Nsim, m);
sim_output.U_corr = zeros(Nsim, m);
sim_output.opt_val = zeros(Nsim, 1);
sim_output.supervision_active = zeros(Nsim, 1);
sim_output.active_region = zeros(Nsim, 1);
sim_output.Nsim = Nsim;

% Construct the high-dimensional dynamical system wrt original space:
Ti = [1 zeros(1,L-1)];
T = [];
Pi = [zeros(L-1,1) eye(L-1); 1 zeros(1,L-1)];
P = [];
for i = 1:m
    T = blkdiag(T,Ti);
    P = blkdiag(P,Pi);
end
[~,~,~,~,~,Pmat,~,~,alpha,Bm,~] = convert2Bru(Ad, Bd, zeros(n,1), zeros(1,n), 0, [], []);
% A_hd_OG = [Ax Au Av]:
A_hd_OG = [ Ad                      Bd                      zeros(n,L*m); 
           -inv(Bm)*alpha*Pmat*Ad  -inv(Bm)*alpha*Pmat*Bd   T; 
            zeros(L*m,n)            zeros(L*m,m)            P];
Ax = A_hd_OG(:,1:n);
Au = A_hd_OG(:,n+1:n+m);
Av = A_hd_OG(:,n+m+1:end);

% Let us solve the following:
% minimize {over u}     (u-u_des)^T (u-u_des)
% s.t.                  Gcl (Ax x_curr + Au u + Av v) < Fcl
for step = 1 : Nsim
    % Get trajectory:
    trj_pos = trajectory.trj_eval(step * dt, 0);
    trj_vel = trajectory.trj_eval(step * dt, 1);
    trj_acc = trajectory.trj_eval(step * dt, 2);
    trj_jrk = trajectory.trj_eval(step * dt, 3);
    % Get nominal input:
    x_ref = [trj_pos; trj_vel; trj_acc];            % reference next state
    sim_output.X_sim(step, :) = x_curr;             % true current state
    err = x_ref - x_curr;                           % error
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
    % Check if the nominal next state belongs to the CIS:
    guard = isContainedImplicit(x_next_des, DSets,[]);
    
    % If desired next state belongs to the CIS we are good:
    if (guard == 1)
        x_curr = x_next_des;
        sim_output.U_corr(step, :) = u_des;
    % Else we correct the nominal input:
    else
        sim_output.supervision_active(step) = step;
        % We solve for the union of CISs by checking if there exists a next
        % state in any of the CISs separately:
        u_cand = zeros(size(Bd,2), Nsets);          % candidate input
        f_cand = zeros(1, Nsets);                   % candidate optimal cost
        v_cand = zeros(size(Bd,2)*(L), Nsets);      % candidate virtual inputs
        parfor k = 1:Nsets            
            % Implicit CIS inequalitites:
            Gstate = DSets{k}.Gstate;
            Ginput = DSets{k}.Ginput;
            Gvirtual = DSets{k}.Gvirtual;
            Fcl = DSets{k}.Fcl;
            Gcl = [Gstate Ginput Gvirtual];
            % wrt to (u,v):
            Aineq = Gcl*[Au Av];
            bineq = Fcl - Gcl * Ax * x_curr;
            
            % Original cost: || u - udes ||^2.
            HOG = blkdiag(eye(size(Ginput,2)), zeros(size(Gvirtual,2)) );
            cOG = [-2 * u_des; zeros(size(Gvirtual,2),1)];
            H = HOG;
            c = cOG;
            
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
                v_cand(:,k) = NaN*ones(3*L,1);
            end

        end
        
        % Sanity check:
        if (min(isnan(f_cand))) 
            disp(step);
            disp(Nsim);
            error("All the Optimizations returned NaN");
        end
        
        % Find minimum cost over all CISs:
        [sim_output.opt_val(step, :), next_cis] = min(f_cand);
        sim_output.active_region(step) = next_cis;
        % Select input corresponding to minimum cost:
        u_corrected = u_cand(:, next_cis);
        sim_output.U_corr(step, :) = u_corrected;
        % Make sure that next state is indeed in the next cis:
        x_next = Ad * x_curr + Bd * u_corrected;
        guard = isContainedImplicit(x_next, DSets, next_cis);
        if (guard==0)
            disp(step);
            disp(Nsim);
            u = u_cand(:,next_cis);
            v = v_cand(:,next_cis);
            Gx = DSets{next_cis}.Gstate;
            Gi = DSets{next_cis}.Ginput;
            Gv = DSets{next_cis}.Gvirtual;
            Gcl = [Gx Gi Gv];
            Fcl = DSets{next_cis}.Fcl;
            sat = (Gcl * [x_curr; u; v] <= Fcl);
            idcs = find(sat==0);
            diff = abs(Gcl(idcs,:) * [x_curr; u; v] - Fcl(idcs));
            disp(max(diff));
            error("Next state not in a CIS!");
        end
        % Update state:
        x_curr = x_next;
    end
    
end