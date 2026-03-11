function [rhsu, rhsv] = rhs_2d_condiff_var(parm, flow)

    % allocation for speed up
    
    rhsu = zeros(parm.m, parm.n);
    rhsv = zeros(parm.m, parm.n);
    
    % computation of right hand side for both equations

    for i = 1 : parm.m
        for j = 1 : parm.n
            
            ad = advection_u_lin2_var(parm, flow, i, j);
            diff = diffusion_u(parm, flow, i, j);
            rhsu(i,j) = -ad + diff;
            
            ad = advection_v_lin2_var(parm, flow, i, j);
            diff = diffusion_v(parm, flow, i, j);
            rhsv(i,j) = -ad + diff;
            
        end
    end

end

% additional functions
%**************************************************************************

function [adv_u] = advection_u_lin2_var(parm, flow, i, j)

    % for convenience

    ip = parm.ip(i);
    im = parm.im(i);
    jp = parm.jp(j);
    jm = parm.jm(j);
    dx = parm.dx;
    dy = parm.dy;
    
    % advection term 
    
    adv_u = (flow.u_tr(i,j) * (flow.u(ip, j) - flow.u(im, j))/(2*dx)) + ...
            (flow.v_tr(i,j) * (flow.u(i, jp) - flow.u(i, jm))/(2*dy));
    

end

function [adv_v] = advection_v_lin2_var(parm, flow, i, j)

    % for convenience

    ip = parm.ip(i);
    im = parm.im(i);
    jp = parm.jp(j);
    jm = parm.jm(j);
    dx = parm.dx;
    dy = parm.dy;
    
    % advection term 
    
    adv_v = (flow.u_tr(i,j) * (flow.v(ip, j) - flow.v(im, j))/(2*dx)) + ...
            (flow.v_tr(i,j) * (flow.v(i, jp) - flow.v(i, jm))/(2*dy));

end




