function [SingleOCP, SequenceOCP] = ConstructNLP(N, x0, Params);
% -------------------------------------------------------------------------
% This file is part of MPC-DICE.
% 
% MPC-DICE -- Model Predictive Control - Dynamic Integrated model of 
%             Climate and Economy
% Copyright (C) 2018 Timm Faulwasser, Christopher M. Kellett, and 
%       Steven R. Weller. 
% 
% MPC-DICE can be downloaded from https://github.com/cmkellett/MPC-DICE.
% 
% MPC-DICE is free software; you can redistribute it and/or modify it 
% under the terms of the GNU Lesser General Public License (LGPL) as 
% published by the Free Software Foundation; either version 3 of the 
% License, or (at your option) any later version.
% 
% MPC-DICE is distributed in the hope that it will be useful, but 
% without any warranty; without even the implied warranty of merchantability 
% or fitness for a particular purpose.  See the GNU Lesser General Public 
% License for more details.                        
%
% -------------------------------------------------------------------------
%
% Construction of NLPs to be solved at each sampling instant
% [SingleOCP, SequenceOCP] = ConstructNLP(N, x0)
% N ... prediction horizon
% x0 ... initial condition
% SingleOCP ... struct with OCP definition for very first time instant 
% SequenceOCP ... struct with OCP definition for MPC solution

% import casadi name space
import casadi.*

%% ========================================================================
% Define data of the DICE OCP
% =========================================================================
nx = 1+6+9+1; % # of states: time, 6 states, 9 auxilliary states, objective
nu = 2; % # of inputs

if length(x0) ~= nx
    warning('Inconsistent initial condition x0')
    warning(['Length(x0) must be ', num2str(nx)])
end


%% ========================================================================
% Construct bounds on states and inputs
% =========================================================================
%u_LB = zeros(nu, N+1);
%u_UB = ones(nu,  N+1);

% To reproduce Nordhaus DICE2013 results, uncomment lines below
u_LB = [zeros(1,N+1); [zeros(1,N-11) Params.optlrsav*ones(1,11) 0]];
u_UB = [[1 1.0*ones(1,27) Params.limmiu*ones(1,N-27)]; ones(1,N+1)];

x_LB = zeros(nx, N+1);
x_LB(4, :) = 1E0; % to prevent numerical problems in logarithmic function 
x_LB(8:14,:) = -inf;
x_LB(nx,:) = -inf; % no lower bound on objective
x_UB = 1E12*ones(nx, N+1); % no upper bound on states
x_UB(end-2:end-1, :) = 1;  % upper bound on shifted states (inputs)
x_UB(2,:) = Params.T_AT_max; % temperature constraint

xu_LB = [x_LB(:); u_LB(:)];
xu_UB = [x_UB(:); u_UB(:)];

% constraints on initial state (first OCP)
x_LB_ini = zeros(nx, N+2);
x_LB_ini(4, :) = 1E0; % to prevent numerical problems in logarithmic function 
x_LB_ini(8:14,:) = -inf;
x_LB_ini(nx,:) = -inf; % no lower bound on objective
x_LB_ini(end-2:end-1, N+2) = [0; 0];
x_UB_ini = 1E12*ones(nx, N+2); % no upper bound on states
x_UB_ini(2,:) = Params.T_AT_max; % temperature constraint

% To reproduce Nordhaus DICE2013 results, uncomment lines below
x_UB_ini(end-2:end-1, :) = Params.limmiu;  % upper bound on shifted states, i.e. inputs
x_UB_ini(end-2:end-1, N+2) = [Params.miu0; 1];  % initial input constraints in
%Nordhaus' code, need to reproduce Nordhaus' results

%x_UB_ini(end-2:end-1, :) = 1.0;  % upper bound on shifted states, i.e. inputs

xu_LB_ini = [x_LB_ini(:); u_LB(:)];
xu_UB_ini = [x_UB_ini(:); u_UB(:)];

%% ========================================================================
% Define NLP
% =========================================================================
% allocate CASADI variables
x           = SX.sym('x', nx, N+1); % states
u           = SX.sym('u', nu, N+1); % inputs
xini        = SX.sym('xini',nx, 1); % initial condition --> parameter of NLP 
eq_con      = SX.sym('eq_con', nx, N+1); % nx * N+1 constraints for the dynamics

x0_ini = set_initial_conditions(1, xini, Params);

% loop over dynamics
j=1;
for  i = 1:N
    if  i == 1
       eq_con(:,1) = x(:,1) - xini; % x(:,1) = x0;       
    end
    
    % equality constraints for dynamics
    eq_con(:,i+1) = x(:,i+1) - dice_dynamics(x(:,i),u(:,i),Params);  
end

% extra constraints, leaving initial states for input shift as decision variables
eq_con_extra = [xini(1:end-3) - x0_ini(1:end-3); xini(end) - x0_ini(end)]; 
eq_con_ini = [eq_con(:); eq_con_extra(:)];

% define the objective (Mayer term)
%obj = ((5 * 0.016408662 * Params.sc.J*x(nx, N+1)) - 3855.106895);
obj = 1*((x(nx, N+1)) ); % --> improved scaling

%% define NLPs 

% NLP for very first OCP
nlp_ini = struct('x', [[x(:); xini(:)]; u(:)], 'f', obj, 'g', eq_con_ini(:));

% NLP for NMPC sequence of OCPs
nlp = struct('x', [x(:); u(:)], 'f', obj, 'g', eq_con(:), 'p', xini);

%% ========================================================================
% construct guess for states and inputs
% =========================================================================

u_guess = .5*ones(nu, N+1);
x_guess = zeros(nx, N+1);
for i = 1:N
    if i == 1
        x_guess(:,1) = x0;
    end
    x_guess(:,i+1) = dice_dynamics(x_guess(:,i),u_guess(:,i),Params);
end
xu_guess = [x_guess(:); u_guess(:)];
xu_guess_ini = [[x_guess(:); x0(:)];u_guess(:)];

%% ========================================================================
% prepare output data
% =========================================================================
SingleOCP.nlp = nlp_ini;
SingleOCP.xu_guess = xu_guess_ini;
SingleOCP.xu_LB = xu_LB_ini;
SingleOCP.xu_UB = xu_UB_ini;

SequenceOCP.nlp = nlp;
SequenceOCP.xu_guess = xu_guess;
SequenceOCP.xu_LB = xu_LB;
SequenceOCP.xu_UB = xu_UB;

end