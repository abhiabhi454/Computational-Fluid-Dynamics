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
% Explicit Euler scheme for time advancement
%
%--------------------------------------------------------------------------
%
% Clear all variables and plots.
format long;
clear; close all;
hold off;


% Set convection velocity
% U0 = 1.0;
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

% Calculate CFL
% session4 slide-12
CFL = U0 * dt / dx;
% disp(['CFL = ', num2str(CFL)]);

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


% Explicit Euler:
%----------------
%
% phinew is phi at new time level
% phinew must be written back to phi for new timestep
%
% Loop over timesteps:
% session4 slide-10
for i = 1 : tsteps

  % Periodic boundary conditions at x=0:
  phinew(1) = phi(1) - U0*dt / (2*dx) * (phi(2) - phi(points-1));


  % Loop over grid points in space:
  for j = 2 : points - 1

    phinew(j) = phi(j) - U0*dt / (2*dx) * (phi(j+1) - phi(j-1));

  end

  % Periodic boundary conditions at x=2*pi:
  phinew(points) = phi(points) - U0*dt / (2*dx) * (phi(2) - phi(points-1));

  % Write new field back to old field:
  phi = phinew;

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
  title(sprintf('explicit Euler | Step: %d | Time: %.2f s', i, i*dt));
  txt = ['U_0 = ', num2str(U0), ', dt = ', num2str(dt), ', dx = ', num2str(dx), ', CFL = ', num2str(CFL, '%.3f')];
  text(0.05, 0.9, txt, 'Units', 'normalized', 'FontSize', 10);
  legend('Numerical (Explicit Euler)', 'Analytical Solution','Location', 'southeast');
  hold off; grid on;
  pause(0.0000001);

% Pause to allow plot capture ()
if mod(i, tsteps/20) == 0
    temp = input('Enter 1 to continue, or 2 to stop: ');
    % ⚠️Clicking the red Stop button (top-right) in MATLAB Online may cause freezing ;(
    if temp == 2
        disp('Simulation stopped');
        break;
    end
end

end

