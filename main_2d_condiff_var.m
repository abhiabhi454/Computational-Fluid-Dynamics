%**************************************************************************
% CFD Lab Summer Semester 2020, Assignment 6
%
% This code solves the following vectorial 2D advection-diffusion problem 
%
% "TWOD_CONDIFF_VAR"
% du/dt = -(u_tr*du/dx + v_tr*du/dy) + nu*(d^2 u/dx^2 + d^2 u/dy^2) 
% dv/dt = -(u_tr*dv/dx + v_tr*dv/dy) + nu*(d^2 v/dx^2 + d^2 v/dy^2)
% where u_tr and v_tr are functions of space
%
% For spatial discretisation a central difference scheme is used
% and for time integration the explicit Euler time integration and a 
% Runge-Kutta method are available.
%
% Periodic boundary condition is applied both in x- and y-direction
%
% authors: D.Quosdorf & Y.Sakai
% june, 2018
%**************************************************************************

% close figures, command window and clear memory
close all; clc; clear

% read infile 
infilename = 'infile_2D_condiff_var.mat';
fprintf('infilename is: %s\n', infilename)

% build structures 'parm' and 'flow'
[parm, flow] = build_structs;
fprintf('struct built\n')

% fill some fields of 'parm' and 'flow' with data from infile
[parm, flow] = set_params(parm, flow, infilename);
fprintf('parameters set\n')

% initialisation of flow field
[parm, flow] = initialize(parm, flow);
fprintf('flow field initialised\n')

% start time integration
[parm, flow] = timeint(parm, flow);
