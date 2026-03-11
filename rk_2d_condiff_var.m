function [flow] = rk_2d_condiff_var(parm, flow)

    % Implement low-storage 3rd-order Runge-Kutta time integration method
    % by filling the indicated missing part

    % 1st Runge-Kutta step
    
    % note: original RHS is computed outside of this function, and stored
    % in flow.rhsu, rhsv
    
    % initialize a 2-D vector S with the size of parm.m and parm.n
    %S = zeros(parm.m, parm.n, 2 );  

    % main R-K
    for i = 1 : parm.m
        for j = 1 : parm.n

            %S(i,j,1) = flow.rhsu(i,j);
            %S(i,j,2) = flow.rhsv(i,j);
            
            %flow.u(i,j) = ???           
            %flow.v(i,j) = ???

            flow.u(i,j) = flow.u(i,j) + parm.dt * flow.rhsu(i,j) / 3;
            flow.v(i,j) = flow.v(i,j) + parm.dt * flow.rhsv(i,j) / 3;
            
        end
    end
    
    % 2nd Runge-Kutta step (note: "secrhs" stands for second RHS)
    
    [flow.secrhsu, flow.secrhsv] = rhs_2d_condiff_var(parm, flow);
    
    for i = 1 : parm.m
        for j = 1 : parm.n
            
            %flow.secrhsu(i,j) = flow.secrhsu(i,j) + ???
            %flow.secrhsv(i,j) = flow.secrhsv(i,j) + ???
            flow.secrhsu(i,j) = flow.secrhsu(i,j) - flow.rhsu(i,j)*5/9;
            flow.secrhsv(i,j) = flow.secrhsv(i,j) - flow.rhsv(i,j)*5/9;
            
            %flow.u(i,j) = flow.u(i,j)+ ???
            %flow.v(i,j) = flow.v(i,j)+ ???
            flow.u(i,j) = flow.u(i,j) + parm.dt * (flow.secrhsu(i,j)) * 15 / 16;
            flow.v(i,j) = flow.v(i,j) + parm.dt * (flow.secrhsv(i,j)) * 15 / 16;
        
        end
    end
    
    % 3rd Runge-Kutta step
    
    [flow.rhsu, flow.rhsv] = rhs_2d_condiff_var(parm, flow);
    
    for i = 1 : parm.m
        for j = 1 : parm.n
            
            %flow.rhsu(i,j) = flow.rhsu(i,j) + ???
            %flow.rhsv(i,j) = flow.rhsv(i,j) + ???
            flow.rhsu(i,j) = flow.rhsu(i,j) - flow.secrhsu(i,j) *153/128;
            flow.rhsv(i,j) = flow.rhsv(i,j) - flow.secrhsv(i,j) * 153/128;
            flow.u(i,j) = flow.u(i,j) + parm.dt * flow.rhsu(i,j) * 8/15;
            flow.v(i,j) = flow.v(i,j) + parm.dt * flow.rhsv(i,j) * 8/15;
            
            %flow.u(i,j) = flow.u(i,j) + ???
            %flow.v(i,j) = flow.v(i,j) + ??? 
        end
    end

end

