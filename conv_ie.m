% Course: CFD Lab
% TU Muenchen, Summer term 2020
%
% This is the Matlab script for the unsteady 1D convection equation
%
% You must fill in the missing parts by yourself!
% Missing parts are marked by ???
%
% Author: Tianshi Sun, tianshi.sun@tum.de
% 
%
%--------------------------------------------------------------------------
%
% time integration of
%
%   d phi            d phi
%  ------- =  -U0 * -------
%    dt               dx
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
% Implicit Euler scheme for time advancement
%
%--------------------------------------------------------------------------
%
% Clear all variables and plots.
format long;
clear;
hold off;

% Set convection velocity
%U0 = 1.0;
U0 = input('Please enter convection velocity U0: ');

% Discrete spacing in space
xend   = 2.0 * pi;
% points = 40; 
points = input('Please enter number of grid points: ');
dx     = xend / ( points - 1 );
% Grid with x locations:
x = 0.0 : dx : xend;

% Discrete spacing in time
% tstep = number of discrete timesteps
tsteps = 1000;
% dt     = 0.1;
dt = input('Please enter time step size dt: ');
tend   = dt * tsteps;

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
% session4 slide-5
a_w = -U0*dt / (2*dx);   
a_p = 1;                     
a_e = U0*dt / (2*dx);  

% Implicit Euler:
%----------------
%
% Loop over timesteps:
for i = 1 : tsteps

  % Periodic boundary conditions at x=0 (left):
  % session4 slide-8,9
  A(1, points-1) = a_w;
  A(1, 1) = a_p;
  A(1, 2) = a_e;
  b(1) = phi(1);

  % Loop over grid points in space:
  % session4 slide-5,6
  for j = 2 : points - 1

    A(j, j-1) = a_w;
    A(j, j) = a_p;
    A(j, j+1) = a_e;
    b(j) = phi(j);

  end

  % Periodic boundary conditions at x=2*pi (right):
  % session4 slide-8,9
  A(points, points-1 ) = a_w;
  A(points, points) = a_p;
  A(points, 2) = a_e;
  b(points) = phi(points);

  % Solve the linear system of equations
  phi = A\b;

  % Analytical solution
  for j = 1 : points
    % ϕ (x, t) = sin (x − U0 t)
    phi_a(j) = sin (x(j) - U0*i*dt );
    % hint: t(i) = i * dt
  end

  % Plot transported wave for each timestep
  plot(x, phi, 'r', x, phi_a, 'g');
  xlabel('x'); xlim([0 2*pi]);
  ylabel('\phi(x, t)');
  title(sprintf('implicit Euler | Step: %d | Time: %.2f s', i, i*dt));
  txt = ['U_0 = ', num2str(U0), ', dt = ', num2str(dt), ', dx = ', num2str(dx)];
  text(0.05, 0.9, txt, 'Units', 'normalized', 'FontSize', 10);
  legend('Numerical (Implicit Euler)', 'Analytical Solution','Location', 'southeast');
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

