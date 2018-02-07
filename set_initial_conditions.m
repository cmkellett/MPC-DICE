function x_zero = set_initial_conditions(varargin)
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
% -------------------------------------------------------------------------
if length(varargin) == 2
     t0 = varargin{1};
     Params = varargin{2};
     
     % initial conditions for auxiliary state variables
     mu0 = Params.miu0; % initial mitigation rate
     s0= 0.259029014481802; %  initial savings rate
     
elseif  length(varargin) == 3 
     t0 = varargin{1};
     x_ini = varargin{2};
     Params = varargin{3};
     
     % initial conditions for auxiliary state variables
     mu0 = x_ini(15);  % initial mitigation rate == casadi var
     s0= x_ini(16); %  initial savings rate == casadi var
end

x0 = [Params.T_AT0 Params.T_LO0 Params.M_AT0 Params.M_UP0 Params.M_LO0 Params.K0]';

Gross_Economic_Output0 = Params.A0*(x0(6)^Params.gamma)*((Params.L0/1000)^(1-Params.gamma));
Emissions0 = 5*(Params.sigma0*(1-mu0)*Gross_Economic_Output0 + Params.EL0);
Damages0 = 1 - (Params.a2*(x0(1)^Params.a3))/(1+Params.a2*(x0(1)^Params.a3));
theta10 = Params.sigma0 * (Params.pb/(1000*Params.theta2)) * (1-Params.deltaPB)^(t0-1);
Net_Economic_Output0 = Damages0 * (1 - theta10*(mu0^Params.theta2))*Gross_Economic_Output0;
C0 = 5*(1-s0) * Net_Economic_Output0;
        
% scaling         
x0(4) = x0(4)/Params.sc.M;
x0(3) = x0(3)/Params.sc.M;
x0(5) = x0(5)/Params.sc.M;
x0(6) = x0(6)/Params.sc.K;
Params.L0 = Params.L0/Params.sc.L;
C0 = C0/Params.sc.C;

x_zero = [t0; x0; Params.sigma0; Params.L0; Params.A0; Params.EL0; Params.f0; ...
    Emissions0; C0; mu0; s0; 0];