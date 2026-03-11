function [flow] = euler_2d_condiff_var(parm, flow)

    for i = 1 : parm.m
        for j = 1 : parm.n
            
            flow.u(i,j) = flow.u(i,j) + parm.dt * flow.rhsu(i,j) * ...
                (1-parm.atbounds(j));
            
            flow.v(i,j) = flow.v(i,j) + parm.dt * flow.rhsv(i,j) * ...
                (1-parm.atbounds(j));
            
        end
    end
       
end

