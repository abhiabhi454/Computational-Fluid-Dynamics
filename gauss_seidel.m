function [flow, num_it] = gauss_seidel(parm, flow, tol)

% compute divergence of intermediate velocity field 
% (note: flow.div is an [m,n]-matrix)
[flow] = div_cal(parm, flow);
fprintf('DEBUG: max div of int vel. is %e\n', div_check(flow));

flow.dp=zeros(parm.m,parm.n);

AP = -parm.dt/parm.rho*(2/(parm.dx)^2+2/(parm.dy)^2);

AN = parm.dt/parm.rho*1/(parm.dy)^2;

AS = parm.dt/parm.rho*1/(parm.dy)^2;

AE = parm.dt/parm.rho*1/(parm.dx)^2;

AW = parm.dt/parm.rho*1/(parm.dx)^2;

err=1;
num_it = 0;

while err>tol && num_it < 1000
    
    temp = flow.dp;
    num_it = num_it + 1;

    for i = 1 : parm.m
            
        ip = parm.ip(i);
        im = parm.im(i);
    
        for j = 1 : parm.n
            
            jp = parm.jp(j);
            jm = parm.jm(j);
            
            flow.dp(i,j)=1/AP*(flow.div(i,j) - ...
                                AS*flow.dp(i,jm) - ...
                                AN*flow.dp(i,jp) - ...
                                AW*flow.dp(im,j) - ...
                                AE*flow.dp(ip,j));
     
            if strcmp(parm.flowtype,'CHANNEL')
                if i==1 && j==1
                   flow.dp(i,j)=0;
                end            
            end

        end

    end
    
    err=max(abs(temp-flow.dp),[],'all');
end

end