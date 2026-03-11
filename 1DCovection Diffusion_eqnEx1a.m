% CFD-Lab: Assignment 1a
% Steady 1D Convection-Diffusion Equation with FD
% TUM - Summer Term 2025
% Group Number: [8]

format long;
clc; close all;

% Parameters
U0 = 1.0;
Gamma = 1.0;
xend = 2.0*pi;

% Discretization
points = 40;
dx = xend / (points - 1);
x = linspace(0, xend, points);

% Boundary Conditions
phi_0 = 0.0;
phi_end = 1.0;

%Initialize
phi_up = zeros(points,1);
phi_cd = zeros(points,1);
phi_analytic = zeros(points,1);
abs_err_cd=zeros(points,1);
abs_err_up=zeros(points,1);
A_up = zeros(points);
A_cd = zeros(points);
b = zeros(points,1);

%Upwind Scheme
for i = 2 : points-1
    a_w = U0/dx + Gamma/dx^2;
    a_p = -U0/dx - 2*Gamma/dx^2;
    a_e = Gamma/dx^2;
    
    A_up(i, i-1) = a_w;
    A_up(i, i)   = a_p;
    A_up(i, i+1) = a_e;
end

% Apply BCs
A_up(1,1) = 1;     b(1) = phi_0;
A_up(points,points) = 1; b(points) = phi_end;

% Solve
phi_up = A_up \ b;
phi_up = phi_up(:);
%Central Difference Scheme
for i = 2 : points-1
    a_w = U0/(2*dx) + Gamma/dx^2;
    a_p = -2*Gamma/dx^2;
    a_e = -U0/(2*dx) + Gamma/dx^2;
    
    A_cd(i, i-1) = a_w;
    A_cd(i, i)   = a_p;
    A_cd(i, i+1) = a_e;
end

A_cd(1,1) = 1;     b(1) = phi_0;
A_cd(points,points) = 1; b(points) = phi_end;

phi_cd = A_cd \ b;
%Analytical Solution
Pe = U0 * x / Gamma;
phi_analytic = (exp(Pe) - 1) / (exp(U0 * xend / Gamma) - 1);
phi_analytic = phi_analytic(:);
phi_cd = phi_cd(:);
err_at_pi_up = abs(phi_up(20)-phi_analytic(20))/phi_analytic(20);
err_at_pi_cd = abs(phi_cd(20)-phi_analytic(20))/phi_analytic(20);
abs_err_up = abs(phi_up - phi_analytic);
abs_err_cd = abs(phi_cd - phi_analytic);
fprintf(['for upwind scheme, since in the incomplete code given on moodle, the number of grid point taken is 40 ' ...
    'therefore x=pi i.e., the midpoint is not included , hence I would ' ...
     'print the relative error at 20th grid point is = : %.2f\n'],err_at_pi_up);
err_at_pi_cd = abs(phi_cd(20)-phi_analytic(20))/phi_analytic(20);
fprintf(['for central difference scheme, since in the incomplete code given on moodle, the number of grid point taken is 40 ' ...
    'therefore x=pi i.e., the midpoint is not included , hence I would ' ...
     'print the relative error at 20th grid point is = : %.2f\n'],err_at_pi_cd);
%Plot Numerical vs Analytical (Upwind)
figure;
plot(x, phi_up, 'r-o', x, phi_analytic, 'k--');
xlabel('x'); ylabel('\phi');
legend('Upwind FD','Analytical');
title('1D Convection-Diffusion: Upwind Scheme');
grid on;

%Plot Numerical vs Analytical (Central Difference)
figure;
plot(x, phi_cd, 'b-o', x, phi_analytic, 'k--');
xlabel('x'); ylabel('\phi');
legend('Central Difference FD','Analytical');
title('1D Convection-Diffusion: Central Difference Scheme');
grid on;

%Absolute Error curve for given grid size of 40
figure;
plot(x, abs_err_up, 'r');hold on
plot(x, abs_err_cd, 'b');
legend('Upwind','Central');
xlabel('x'); ylabel('|\phi_{num} - \phi_{analytic}|');
title('Absolute Error Distribution');
grid on;

%---------------------
%Test Cases (Tasks 2–5)
test_cases = [
    -10, 5;
    10, 5;
    -10, 51;
    10, 51
    ];

figure;
case_num = 1;

for k = 1:size(test_cases,1)
    u0 = test_cases(k,1);
    points = test_cases(k,2);
    dx = xend / (points - 1);
    x = linspace(0, xend, points);
    Gamma = 1.0;
    
    % Analytical solution
    phi_exact = (exp(u0 * x / Gamma) - 1) / (exp(u0 * xend / Gamma) - 1);
    
    % --- Upwind ---
    A = zeros(points); b = zeros(points,1);
    for i = 2 : points-1
        a_w = u0/dx + Gamma/dx^2;
        a_p = -u0/dx - 2*Gamma/dx^2;
        a_e = Gamma/dx^2;
        
        A(i,i-1) = a_w;
        A(i,i)   = a_p;
        A(i,i+1) = a_e;
    end
    A(1,1) = 1; b(1) = 0.0;
    A(points,points) = 1; b(points) = 1.0;
    
    phi_up = A\b;
    
    subplot(4,2,2*k-1)
    plot(x, phi_up, 'r-o', x, phi_exact, 'k--');
    title(sprintf('Case %d Upwind: U0=%.1f, N=%d', case_num, u0, points));
    xlabel('x'); ylabel('\phi'); legend('Upwind', 'Analytical'); grid on;
    
    % --- Central ---
    A = zeros(points); b = zeros(points,1);
    for i = 2 : points-1
        a_w = u0/(2*dx) + Gamma/dx^2;
        a_p = -2*Gamma/dx^2;
        a_e = -u0/(2*dx) + Gamma/dx^2;
        
        A(i,i-1) = a_w;
        A(i,i)   = a_p;
        A(i,i+1) = a_e;
    end
    A(1,1) = 1; b(1) = 0.0;
    A(points,points) = 1; b(points) = 1.0;
    
    phi_cd = A\b;
    
    subplot(4,2,2*k)
    plot(x, phi_cd, 'b-s', x, phi_exact, 'k--');
    title(sprintf('Case %d Central: U0=%.1f, N=%d', case_num, u0, points));
    xlabel('x'); ylabel('\phi'); legend('Central', 'Analytical'); grid on;
    
    case_num = case_num + 1;
end


% -----------------------------
% Convergence Study (Upwind & Central)
% -----------------------------
grid_sizes = [11, 21, 31, 41, 51, 61];  % Ensure x = pi is a grid point
errors_up = zeros(size(grid_sizes));
errors_cd = zeros(size(grid_sizes));
dx_values = zeros(size(grid_sizes));
rms_cd = zeros(size(grid_sizes));
rms_up = zeros(size(grid_sizes));
for k = 1:length(grid_sizes)
    points = grid_sizes(k);
    dx = xend / (points - 1);
    x = linspace(0, xend, points);
    %[~, idx_pi] = min(abs(x - pi));  % Ensure x = pi is captured
    idx_pi = ((grid_sizes(k)-1)/2)+1;
    % Analytical solution
    phi_exact = (exp(U0 * x / Gamma) - 1) / (exp(U0 * xend / Gamma) - 1);
    phi_exact = phi_exact(:);
    %% --- Upwind ---
    A = zeros(points, points);
    b = zeros(points, 1);
    for i = 2:points-1
        a_w = U0/dx + Gamma/dx^2;
        a_p = -U0/dx - 2*Gamma/dx^2;
        a_e = Gamma/dx^2;
        A(i, i-1) = a_w;
        A(i, i)   = a_p;
        A(i, i+1) = a_e;
    end
    A(1,1) = 1; b(1) = phi_0;
    A(end,end) = 1; b(end) = phi_end;
    phi_num_up = A \ b;
    phi_num_up = phi_num_up(:);
    errors_up(k) = (abs(phi_exact(idx_pi) - phi_num_up(idx_pi))/phi_exact(idx_pi));
    %error_abs_up = abs(phi_exact - phi_num_up);
    %size(error_abs_up)
    %rms_up(k) = sqrt(mean(error_abs_up.^2));
    %% --- Central Difference ---
    A = zeros(points, points);
    b = zeros(points, 1);
    for i = 2:points-1
        a_w = U0/(2*dx) + Gamma/dx^2;
        a_p = -2*Gamma/dx^2;
        a_e = -U0/(2*dx) + Gamma/dx^2;
        A(i, i-1) = a_w;
        A(i, i)   = a_p;
        A(i, i+1) = a_e;
    end
    A(1,1) = 1; b(1) = phi_0;
    A(end,end) = 1; b(end) = phi_end;
    phi_num_cd = A \ b;
    phi_num_cd = phi_num_cd(:);
    errors_cd(k) = (abs(phi_exact(idx_pi) - phi_num_cd(idx_pi))/phi_exact(idx_pi));
    %error_abs_cd = abs(phi_exact - phi_num_cd);
    %rms_cd(k) = sqrt(mean(error_abs_cd.^2));
    dx_values(k) = dx;
end

figure;
plot(log(dx_values),log(errors_up),'o-','LineWidth', 2);hold on;
plot(log(dx_values),log(errors_cd),'s-','LineWidth', 2);
legend('Upwind', 'Central Difference');
grid on;
xlabel('log(dx)'); ylabel('log(Relative Error at x = \pi)');
legend('Upwind', 'Central Difference');
title('Convergence Study: Error at x= \pi vs grid spacing h');

% Estimate and print order of accuracy
p_up = polyfit((dx_values),(errors_up), 1);
p_cd = polyfit((dx_values),(errors_cd), 1);
fprintf('Estimated order of accuracy:\n');
fprintf('  Upwind: %.2f (expected ~1)\n', p_up(1));
fprintf('  Central: %.2f (expected ~2)\n', p_cd(1));

