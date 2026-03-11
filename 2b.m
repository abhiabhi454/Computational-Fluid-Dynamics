% Course: CFD Lab
% TU-Muenchen, SS 2020
%
% This is the Matlab script for the unsteady 1D convection diffusion equation
%
% You must fill in the missing parts by yourself!
% Missing parts are marked by ???
%
% Tianshi Sun, tianshi.sun@tum.de
% 
%
%--------------------------------------------------------------------------
%
% time integration of
%
%   d phi            d phi             d^2 phi
%  ------- =  -U0 * ------- + Gamma * ---------
%    dt               dx                dx^2 
%
% 0 <= x <= 2pi
%
% periodic boundary condition
% phi(0) = phi(2pi)
%
% initial condition
% t = 0  ==>  phi = sin(x)
%
% Central difference scheme (CDS) for spatial discretization
% Crank Nicolson scheme for time advancement
%
%--------------------------------------------------------------------------
%
% Clear all variables and plots.
format long;
clear;
hold off;

% Set convection velocity
% U0 = 1.0;
U0 = input('Please enter convection velocity U0: ');

% Set diffusion coefficient
% Gamma = 1.0;
Gamma = input('Please enter diffusion coefficients Γ: ');

% Discrete spacing in space
xend   = 2.0 * pi;
%points = 40; 
points = input('Please enter number of grid points: ');
dx     = xend / ( points - 1 );
% Grid with x locations:
x = 0.0 : dx : xend;

% Discrete spacing in time
% tstep = number of discrete timesteps
tsteps = 1000;
%dt     = 0.1;
dt = input('Please enter time step size dt: ');
tend   = dt * tsteps;

% Diffusion Number
% Session5 slides-9
D = Gamma * dt / dx^2;

% Initialise coefficient matrix A, constant vector b
% and solution vector phi
A   = zeros(points,points);
b   = zeros(points,1);
phi = zeros(points,1);

% Initialise the solution (initial condition)
% Loop over grid points in space:
for j = 1 : points
   phi(j) = sin(x(j));
end

% Check initial field:
plot(x, phi, 'r');
xlabel('x'); xlim([0 2*pi]);
ylabel('\phi(x, 0)');
title('Initial Condition: \phi(x, 0) = sin(x)');
legend('Initial Field');
hold on;
pause(3);

% Compute coefficients of matrix A
% session5 slide-5
a_w = U0 / (2*dx) + Gamma / dx^2;
a_p = -2*Gamma / dx^2;
a_e = -U0 / (2*dx) + Gamma / dx^2;

% Crank Nicolson:
%----------------
%

% Coefficient A at n+1
% Session5 slide-6
AW = -dt/2 * a_w;
AP = 1 - dt/2 * a_p;
AE = -dt/2 * a_e;

% Loop over timesteps:
for i = 1 : tsteps

  % Periodic boundary conditions at x=0:
  % session5 slide-6~8
  A(1,points-1) = AW;
  A(1,1) = AP;
  A(1,2) = AE;
  b(1) = dt/2 * (a_w*phi(points - 1) + a_p*phi(1) + a_e*phi(2)) + phi(1);

  % Loop over grid points in space:
  for j = 2 : points - 1

    A(j, j-1) = AW;
    A(j, j)   = AP;
    A(j, j+1) = AE;
    b(j) = dt/2 * (a_w * phi(j-1) + a_p * phi(j) + a_e * phi(j+1)) + phi(j);

  end

  % Periodic boundary conditions at x=2*pi:
  A(points, points-1) = AW;
  A(points, points)   = AP;
  A(points, 2) = AE;
  b(points) = dt/2 * (a_w * phi(points-1) + a_p * phi(points) + a_e * phi(2)) + phi(points);

  % Solve the linear system of equations
  phi = A\b;

  % Analytical solution
  % ϕ(x, t) = exp(-Gamma t) * sin(x - U0 t)
  for j = 1 : points
    phi_a(j) = exp(-Gamma * i * dt) * sin(x(j) - U0*i*dt);
    % hint: t(i) = i * dt
  end

  % Plot transported wave for each timestep
  plot(x, phi, 'r', x, phi_a, 'g');
  xlabel('x'); xlim([0 2*pi]);
  ylabel('\phi(x, t)');
  title(sprintf('Crank-Nicolson | Step: %d | Time: %.2f s', i, i*dt));
  txt = ['U_0 = ', num2str(U0), ', Γ = ', num2str(Gamma) ', dt = ', num2str(dt), ', dx = ', num2str(dx), ', Diffusion Number = ', num2str(D, '%.3f')];
  text(0.05, 0.9, txt, 'Units', 'normalized', 'FontSize', 10);
  legend('Numerical (Explicit Euler)', 'Analytical Solution','Location', 'southeast');
  hold off; grid on;
  pause(0.0000001);

% Pause(every 50 steps) to allow plot capture
if mod(i, tsteps/20) == 0
    temp = input('Enter 1 to continue, or 2 to stop: ');
    % ⚠️Clicking the red Stop button (top-right) in MATLAB Online may cause freezing ;(
    if temp == 2
        disp('Simulation stopped');
        break;
    end
end
end
