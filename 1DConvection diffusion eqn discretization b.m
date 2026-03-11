% Course: Computational Fluid Dynamics 
% Semester: SoSe2020
%% Group Number: [8]
% Assignment 3 
%
% This is a Matlab code to solve 1D steady advection-diffusion equation
% discritised by finite-volume schemes 
%
% differential form of advection-diffusion eqn.:
%
%       d phi          d^2 phi
%  -U0 ------- + Gamma ------- = 0
%        dx             dx^2
%
% You are asked to fill in the missing parts to complete the implementation.
% Missing parts are marked by ???
%
% Author: Yoshiyuki Sakai
% Email: yoshiyuki.sakai@tum.de


% Clear all variables and plots.
format long;
clear;
hold off;

% Set constant advection velocity
U0 = 1.0;

% Set constant diffusivity
Gamma = 1.0;

% Set up grid cells
xend = 2.0 * pi;
cells = 51; 
dx = xend /(cells - 1);


% Array of grid cell centre locations:
x = linspace(0, xend, cells);

% Initialization of cell-averaged field
phi = zeros(cells,1);

% Initialization of matrix A
A = zeros(cells,cells);

% Initialization of vector b
b = zeros(cells,1);

% Boundary cell face values
phi_0   = 0.0;
phi_end = 1.0;

% Loop over grid cells
% NB: boundary cells are excluded
for i = 2 : cells-1

     a_w = (U0/2) + (Gamma/dx);
     a_p = -(2*Gamma)/dx;
     a_e = -(U0/2)+ (Gamma/dx);
     
%     assign values to LHS matrix A
     A(i,i) = a_p;
     A(i,i-1) = a_w;
     A(i,i+1) = a_e;

%     assign values to RHS vector b
     b(i) = phi_0;
end

% Boundary conditions (Dirichlet at boundary cell faces)

% at i = 1
 A(1,1) = -(U0/2)-3*(Gamma/dx);
 A(1,2) = a_e;
 b(1) = -(U0+2*(Gamma/dx))*phi_0;

% at i = cells
 A(cells,cells) = -(3*Gamma)/dx;
 A(cells,cells-1) = (Gamma/dx) + U0/2;
 b(cells) = ((U0/2) - (2*Gamma)/dx)*phi_end;

% Solution of the linear system
phi = A\b;
phi = phi(:);
% Compute analytical solution
Pe = U0 * x / Gamma;
phi_analytic = (exp(Pe) - 1) / (exp(U0 * xend / Gamma) - 1);
phi_analytic = phi_analytic(:);
% Compute relative error
nn = ceil(cells / 2);
err_rel = abs(phi_analytic(nn) - phi(nn)) / phi_analytic(nn);
fprintf(['The relative error at x=pi is equal to: %.4f\n'],err_rel);
% Compute mean error
phi_diff = phi_analytic - phi;
phi_diff = phi_diff(:);
err_mean = sqrt(mean(phi_diff.^2)) / mean(phi_analytic);
% Plot the numerical and analytical solutions
%figure;
plot(x, phi, 'b-o', x, phi_analytic, 'k--');
xlabel('x'); ylabel('\phi');
legend('Numerical','Analytical');
title('1D Convection-Diffusion: Numerical Vs Analytical');
grid on;

%% -----Task 2: Test Cases from experiment 1a --------

%Test Cases
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
    cells = test_cases(k,2);
    dx = xend / (cells - 1);
    x = linspace(0, xend, cells);
    Gamma = 1.0;
    
    % Analytical solution
    phi_exact = (exp(u0 * x / Gamma) - 1) / (exp(u0 * xend / Gamma) - 1);
    
    
    % --- Numerical ---
    A = zeros(cells); b = zeros(cells,1);
    for i = 2 : cells-1
        a_w = (u0/2) + (Gamma/dx);
        a_p = -(2*Gamma)/dx;
        a_e = -(u0/2)+ (Gamma/dx);
        
        A(i,i-1) = a_w;
        A(i,i)   = a_p;
        A(i,i+1) = a_e;
    end


    % Boundary conditions (Dirichlet at boundary cell faces)
    % at i = 1

    A(1,1) = -(u0/2)-3*(Gamma/dx);
    A(1,2) = a_e;
    b(1) = -(u0+2*(Gamma/dx))*phi_0;

    % at i = cells
    A(cells,cells) = -(3*Gamma)/dx;
    A(cells,cells-1) = (Gamma/dx) + u0/2;
    b(cells) = ((u0/2) - (2*Gamma)/dx)*phi_end;
    
    
    phi_num = A\b;
    
    subplot(2,2,k)
    plot(x, phi_num, 'b-s', x, phi_exact, 'k--');
    title(sprintf('Case %d Numerical Solution: U0=%.1f, N=%d', case_num, u0, cells));
    xlabel('x'); ylabel('\phi'); legend('Central', 'Analytical'); grid on;
    
    case_num = case_num + 1;
    end


%% ------Task 3 and 4-----With differnet Grid size------------

% Grid resolution test

cell_sizes = [25, 51, 101];
dx_list = zeros(size(cell_sizes));
err_rel = zeros(size(cell_sizes));
err_mean = zeros(size(cell_sizes));
for idx = 1:length(cell_sizes)

    cells = cell_sizes(idx);
    dx = xend / (cells - 1);  % Number of *nodes* is cells; intervals = cells - 1
    dx_list(idx) = dx;

    x = linspace(0, xend, cells)';  % Column vector

    % Initialize matrices
    A = zeros(cells);
    b = zeros(cells, 1);
    phi = zeros(cells, 1);

    % Boundary conditions
    phi_0 = 0.0;
    phi_end = 1.0;

    % Loop over interior cells (central differencing)
    for i = 2 : cells - 1
        a_w = (U0/2) + (Gamma/dx);
        a_p = -(2*Gamma)/dx;
        a_e = -(U0/2)+ (Gamma/dx);

        A(i, i-1) = a_w;
        A(i, i)   = a_p;
        A(i, i+1) = a_e;

        b(i) = phi_0;  % RHS is zero inside
    end

    % Boundary at x = 0 (Dirichlet)
    
    % at i = 1
    A(1,1) = -(U0/2)-3*(Gamma/dx);
    A(1,2) = a_e;
    b(1) = -(U0+2*(Gamma/dx))*phi_0;

    % Boundary at x = 2*pi (Dirichlet)
   
    A(cells,cells) = -(3*Gamma)/dx;
    A(cells,cells-1) = (Gamma/dx) + U0/2;
    b(cells) = ((U0/2) - (2*Gamma)/dx)*phi_end;

    % Solve linear system
    phi = A \ b;
    phi = phi(:);
    % Analytical solution
    phi_exact = (exp(U0 * x / Gamma) - 1) / (exp(U0 * xend / Gamma) - 1);
    phi_exact = phi_exact(:);
    % Relative error at x = pi
    idx_pi = ceil(cells / 2);
    Erel(idx) = abs(phi_exact(idx_pi) - phi(idx_pi)) / phi_exact(idx_pi);
    % Mean error (normalized RMS)
    phi_diff = phi_exact - phi;
    phi_diff = phi_diff(:);
    Emean(idx) = sqrt(mean(phi_diff.^2)) / mean(phi_exact);
end


% Plot the error as function of dx in log-log scale

% Plot Erel vs dx (log-log)
figure;
loglog(dx_list, Erel, 'ro-', 'LineWidth', 2);
xlabel('Grid spacing (dx)');
ylabel('Relative Error at x = \pi');
title('E_{rel} vs dx');
grid on;

% Plot Emean vs dx (log-log)
figure;
loglog(dx_list, Emean, 'bs-', 'LineWidth', 2);
xlabel('Grid spacing (dx)');
ylabel('Normalized RMS Error (E_{mean})');
title('E_{mean} vs dx');
grid on;

% Estimate order of accuracy (slope)
p_rel = polyfit(log(dx_list), log(Erel), 1);
p_mean = polyfit(log(dx_list), log(Emean), 1);

fprintf('Order of accuracy (E_rel):  %.2f (expected ~2 for central scheme)\n', p_rel(1));
fprintf('Order of accuracy (E_mean): %.2f (expected ~2 for central scheme)\n', p_mean(1));
