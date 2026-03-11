function [parm, flow, nums_it, time_per_timestep] = timeint(parm, flow, solver, tol)
    if nargin < 4
        tol = 1e-4; % default tolerance if not provided
    end

    nums_it = zeros(1, parm.ntst);
    time_per_timestep = zeros(1, parm.ntst); % array to store timing per timestep

    % loop over timesteps 
    for itst = 1 : parm.ntst
        
        tic; % start timing
        
        disp(['time step ', num2str(itst)]);

        % compute right hand side of Navier-Stokes
        % (note: pressure gradient term is practically neglected)
        [flow.rhsu, flow.rhsv] = rhs_ns(parm, flow);

        % compute provisional velocity field 
        % (note: flow.u and flow.v are updated)
        [flow] = runge_kutta_2d_vec(parm, flow);

        % compute pressure difference (dp) field based on the
        
        switch solver

            case {"direct"}
                % provisional velocity (direct pressure correction)
                [flow] = direct_press_corr(parm, flow);

            case {"iterative"}
                % Gauss Seidel (iterative) pressure correction
                [flow, num_it] = gauss_seidel(parm, flow, tol);
                nums_it(itst) = num_it;
                fprintf("Timestep %d required %d iterations\n", itst, num_it);
        end

        % finally correct the velocity field by projection
        % (note: pressure field is also updated)
        [flow] = project(parm, flow);

        % check if the updated velocity fields are indeed divergence-free
        [flow] = div_cal(parm, flow);
        fprintf('DEBUG: max div of final vel. is %e\n', div_check(flow));

        if strcmp(parm.flowtype, 'CHANNEL')
            flow.u_c(itst) = flow.u(ceil(parm.m/2), ceil(parm.n/2));
            flow.u_c_profile(itst, :) = flow.u(ceil(parm.m/2), :);
        end
        flow.u_t(itst, :, :) = transpose(flow.u);
        flow.v_t(itst, :, :) = transpose(flow.v);

        time_per_timestep(itst) = toc; % record elapsed time for this timestep
    end    
end
