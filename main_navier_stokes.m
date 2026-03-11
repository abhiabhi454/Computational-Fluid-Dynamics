%**************************************************************************
% CFD Lab Summer Semester 2025, Assignment 5
% This code solves the 2D unsteady Navier-Stokes equation 
% Central difference (space) + Runge-Kutta (time)
% Gauss-Seidel & Direct pressure correction
% Group 8
%**************************************************************************

close all; clc; clear

resolutions = [20, 40, 60, 80, 100];
solvers = ["direct", "iterative"];
Ntst = 1000;

avg_timestep_time = zeros(2, length(resolutions));  % [solver, resolution]

% Toggle modes
running_cost_comparison = false;
tolerance_comparison = false;

if running_cost_comparison
    for ii = 1:length(resolutions)
        for jj = 1:2
            infilename = 'infile_channel.mat';
            fprintf('infilename is: %s\n', infilename);

            [parm, flow] = build_structs;
            [parm, flow] = set_params(parm, flow, infilename);

            parm.m = resolutions(ii) + 1;
            parm.n = resolutions(ii) + 1;
            parm.dt = 0.0025 * (20 / resolutions(ii))^2;
            parm.ntst = Ntst;

            [parm, flow] = initialize(parm, flow, solvers(jj));
            fprintf('flow field initialised\n');

            % Measure time per timestep using updated timeint
            [parm, flow, ~, time_per_timestep] = timeint(parm, flow, solvers(jj));

            avg_timestep_time(jj, ii) = mean(time_per_timestep);
        end
    end

    % Plot cost comparison
    figure;
    hold on;
    loglog(resolutions, avg_timestep_time(1,:), 'r-o', 'LineWidth', 1.5);
    loglog(resolutions, avg_timestep_time(2,:), 'b-o', 'LineWidth', 1.5);
    grid on;
    legend("Direct Solver", "Gauss-Seidel (Iterative)", 'Location', "northwest");
    xlabel("Grid size m = n");
    ylabel("Average time per timestep (s)");
    title('Solver Time per Timestep vs Grid Size');
    set(gca, 'FontSize', 12);
    hold off;

elseif tolerance_comparison
    solver = "iterative";
    infilename = 'infile_channel.mat';

    [parm, flow] = build_structs;
    [parm, flow] = set_params(parm, flow, infilename);

    tols = [1e-4, 1e-5, 1e-6];
    iterations = zeros(3, parm.ntst);

    for tt = 1:3
        [parm, flow] = initialize(parm, flow, solver);
        [parm, flow, nums_it] = timeint(parm, flow, solver, tols(tt));
        iterations(tt,:) = nums_it;
    end

    figure("Position", [0,0,600,400])
    semilogx(1:parm.ntst, iterations(1,:), 'r', 'LineWidth', 1.5)
    hold on
    semilogx(1:parm.ntst, iterations(2,:), 'g', 'LineWidth', 1.5)
    semilogx(1:parm.ntst, iterations(3,:), 'b', 'LineWidth', 1.5)
    grid on

    ylim([0 250])
    xticks(10.^(1:5))
    xlim([10^1, 10^5])

    legend("tol = 1e-4", "tol = 1e-5", "tol = 1e-6", 'Location', "northeast")
    xlabel("Time step (log scale)")
    ylabel("Iterations per time step")
    title("Gauss-Seidel Iterations vs Time Step")
    set(gca, 'FontSize', 12)
    hold off

else
    % Single run with active plotting
    solver = "iterative";  % or "direct"
    infilename = 'infile_channel.mat';

    [parm, flow] = build_structs;
    [parm, flow] = set_params(parm, flow, infilename);

    parm.m = 21;
    parm.n = 11;
    parm.dt = 0.0025;
    parm.nu = 0.5;
    parm.rho = 1;
    parm.ntst = 200;

    [parm, flow] = initialize(parm, flow, solver);
    [parm, flow, nums_it] = timeint(parm, flow, solver, 1e-4);

    mean_its = mean(nums_it);
    fprintf("Mean Gauss-Seidel iterations per timestep: %.2f\n", mean_its);

    % Plot flow result
    x = 0:parm.dx:parm.xl;
    y = 0:parm.dy:parm.yl;
    figure;
    contourf(x, y, flow.p');
    hold on;
    quiver(x, y, flow.u', flow.v');
    colorbar;
    title("Final Pressure and Velocity Field");
    xlabel("x");
    ylabel("y");
end
