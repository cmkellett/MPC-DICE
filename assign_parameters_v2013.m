function Params = assign_parameters_v2013
% -------------------------------------------------------------------------
%  This file is part of MPC-DICE.
% 
% MPC-DICE -- Model Predictive Control - Dynamic Integrated model of 
%             Climate and Economy
% Copyright (C) 2018 Timm Faulwasser, Christopher M. Kellett, and 
%       Steven R. Weller. 
% 
% MPC-DICE can be downloaded from https://github.com/cmkellett/DICE2013R-mc.
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
% 
% -------------------------------------------------------------------------

%% ------------------------------------------------------------------------
% The following is the unique identifier that will be used to label 
% the saved file from the model run
Params.parameter_set = 'v2013';

% Limit on maximum atmospheric temperature
Params.T_AT_max = 20; 

% -------------------------------------------------------------------------
% Constants
    Params.N = 60;                  % Horizon length
    Params.BaseYear = 2010;        % First calendar year
    Params.eta = 3.8;              % Forcings of equilibrium CO2 doubling (GAMS fco22x)
    Params.M_AT_Base = 588;        % Base atm carbon concentration (GAMS in FORC(t) eqn)
    Params.deltaK = 0.1;            % Capital depreciation (5 year) (GAMS dk)
    Params.gamma = 0.3;            % Capital elasticity in production function (GAMS gama)
    Params.theta2 = 2.8;           % Exponent of control cost function (GAMS expcost2)
    Params.a2 = 0.00267;           % Damage multiplier
    Params.a3 = 2;                 % Damage exponent
    Params.alpha = 1.45;           % Elasticity of marginal utility of consumption (GAMS elasmu)
    Params.rho = 0.015;            % Initial rate of social time preference per year (GAMS prstp)
    Params.xi1 = 0.098;            % Climate equation coefficient for upper level (GAMS c1)
    Params.xi2 = 3/11;             % Conversion factor from GtC to CtCO2

    Params.limmiu = 1.2;           % Upper limit on emissions drawdown

    % Climate Model Diffusion Parameters
    %   Note Vanilla DICE2013R uses nonstandard column-row indexing, where
    %   DICE2013R-mc uses standard row-column indexing.  The diffusion
    %   parameters are calculated based on other parameters and we include
    %   those here as they change between 2013 and 2016 versions.
    c3 = 0.088;
    c4 = 0.025;
    t2xco2 = 2.9;

    phi11 = 1-Params.xi1*((Params.eta/t2xco2) + c3);
    phi12 = Params.xi1*c3;
    phi21 = c4;
    phi22 = 1-c4;

    Params.Phi_T = [phi11 phi12; phi21 phi22];

    % Carbon Cycle Model Diffusion Parameters 
    %   Note Vanilla DICE2013R uses nonstandard column-row indexing, where
    %   DICE2013R-mc uses standard row-column indexing.  The diffusion
    %   parameters are calculated based on other parameters and we include
    %   those here as they change between 2013 and 2016 versions.
    b12 = 0.088;
    b23 = 0.0025;
    mateq = 588;
    mueq = 1350;
    mleq = 10000;

    zeta11 = 1 - b12;
    zeta21 = b12;
    zeta12 = (mateq/mueq)*zeta21;
    zeta22 = 1 - zeta12 - b23;
    zeta32 = b23;
    zeta23 = zeta32*(mueq/mleq);
    zeta33 = 1 - zeta23;

    Params.Phi_M = [zeta11 zeta12 0; zeta21 zeta22 zeta23; 0 zeta32 zeta33];

    % -------------------------------------------------------------------------
    % Exogenous Signal Constants
    Params.L0 = 6838;           % Initial population
    Params.La = 10500;          % Asymptotic population
    Params.lg = 0.134;          % Population growth rate

    Params.EL0 = 3.3;           % Initial land use emissions
    Params.deltaEL = 0.2;       % Land use emissions decrease rate

    Params.A0 = 3.8;            % Initial Total Factor Productivity (TFP)
    Params.ga = 0.079;          % Initial TFP rate
    Params.deltaA = 0.006;      % TFP increase rate

    Params.pb = 344;            % Initial backstop price
    Params.deltaPB = 0.025;     % Decline rate of backstop price

    e0 = 33.61;                 % Initial emissions
    q0 = 63.69;                 % Initial global output
    Params.miu0 = 0.039;               % Initial mitigation rate
    Params.sigma0 = e0/(q0*(1-Params.miu0)); % Calculated initial emissions intensity
    Params.gsigma = 0.01;      % Emissions intensity base rate
    Params.deltasigma = 0.001; % Decline rate of emissions intensity

    Params.f0 = 0.25;           % Initial forcings of non-CO2 GHGs
    Params.f1 = 0.7;            % Forcings of non-CO2 GHGs in 2100
    Params.tforce = 18;         % Slope of non-CO2 GHG forcings

    Params.optlrsav = (Params.deltaK + .004)/(Params.deltaK ...
        + 0.004*Params.alpha + Params.rho) * Params.gamma;
    
    % initial condition for states
    Params.T_AT0 = 0.8;
    Params.T_LO0 = 0.0068;
    Params.M_AT0 = 830.4;
    Params.M_UP0 = 1527;
    Params.M_LO0 = 10010;
    Params.K0 = 135;
    
%% -----------------------------------------------------------------------
%% DO NOT EDIT BELOW THIS LINE UNLESS YOU REALLY KNOW WHAT YOU'RE DOING!!!
%% -----------------------------------------------------------------------

%% scaling of selected state variables
% required for numerical stability
sc.L        = 1E2;
sc.M        = 1E2;
sc.J        = 1E4;
sc.K        = 1E2;
sc.C        = 1E2;
Params.sc   = sc;

%% Options for IPOPT
opts = struct;
opts.ipopt.max_iter    = 3000;
opts.ipopt.print_level = 5;%0,3
opts.print_time        = 0;
opts.ipopt.acceptable_tol =1e-8;
opts.ipopt.acceptable_obj_change_tol = 1e-8;

%opts.ipopt.acceptable_constr_viol_tol = 1E-6;
%opts.ipopt.compl_inf_tol = 1E-6;
%opts.ipopt.dual_inf_tol    = 1E-1;  
%opts.ipopt.mu_strategy                 = 'adaptive';
%opts.ipopt.obj_constr_filter = 'never-monotone-mode'%'kkt-error';
%opts.ipopt.least_square_init_duals     = 'yes';
%opts.ipopt.least_square_init_primals     = 'yes';
%opts.ipopt.bound_relax_factor = 0;
%opts.ipopt.alpha_for_y = 'safer-min-dual-infeas';
%opts.ipopt.alpha_for_y = 'primal-and-full';
%opts.ipopt.alpha_for_y_tol = 5;

%opts.ipopt.hessian_approximation = 'limited-memory';

Params.opts = opts;

