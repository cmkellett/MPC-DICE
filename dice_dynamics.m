function f = dice_dynamics(x,u,Params)
% -------------------------------------------------------------------------
%  This file is part of MPC-DICE.
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
%
% -------------------------------------------------------------------------
% f = dice_dynamics(x,u,Params)
% x ... current state
% u ... current input
% Params .. parameters struct
%
% -------------------------------------------------------------------------
%                               Definition of State Variables
% -------------------------------------------------------------------------
% i = x(1); % time index considered as state variable
% T_AT = x(2); 
% T_LO = x(3);
% M_AT = x(4);
% M_UP = x(5);
% M_LO = x(6);
% K = x(7);
%
% Time-varying parameters reformulated as additional states
% sigma = x(8);
% L = x(9);
% A_TFP = x(10);
% E_LAND = x(11);
% F_EX = x(12);
% E = x(13);  emissions
% C = x(14);  consumption
% mu = x(15); mitigation rate at time i
% s = x(16);  savings rate at time i
% 
% J = x(17); objective
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
%                               Definition of Control Variables
% -------------------------------------------------------------------------
% mu_NEXT = u(1); mitigation rate at time i+1
% s_NEXT = u(2);  savings rate at time i+1
% -------------------------------------------------------------------------


%% ------------------------------------------------------------------------
%                               Definition of State Variables
% -------------------------------------------------------------------------
sc = Params.sc; % scaling parameters for states

i    = x(1); % time index considered as state variable
T_AT = x(2); 
T_LO = x(3);
M_AT = sc.M*x(4);
M_UP = sc.M*x(5);
M_LO = sc.M*x(6);
K    = sc.K*x(7);

% Time-varying parameters reformulated as additional states
sigma = x(8);
L     = sc.L*x(9);
A_TFP = x(10);
E_LAND = x(11);
F_EX   = x(12);
E     = x(13);
C     = sc.C*x(14);

% Controls reformulated as states
mu = x(15);
s  = x(16);
 
% Objective
J = sc.J*x(17);  % scale obj 

T = [T_AT; T_LO];
M = [M_AT; M_UP; M_LO];

%% ------------------------------------------------------------------------
%                               Controls
% -------------------------------------------------------------------------

% Shifted Mitigation Rate
mu_NEXT = u(1);

% Shifted Savings Rate
s_NEXT = u(2);

%% ------------------------------------------------------------------------
%                           Unpack Parameters
% -------------------------------------------------------------------------
eta         = Params.eta;
M_AT_Base   = Params.M_AT_Base;
deltaK      = Params.deltaK;
gamma       = Params.gamma;
theta2      = Params.theta2;
alpha       = Params.alpha;
rho         = Params.rho;
xi1         = Params.xi1;
xi2         = Params.xi2;
Phi_T       = Params.Phi_T;
Phi_M       = Params.Phi_M;
zeta11      = Phi_M(1,1);
zeta12      = Phi_M(1,2);
a2          = Params.a2;
a3          = Params.a3;

La          = Params.La;          % Asymptotic population
lg          = Params.lg;          % Population growth rate

EL0         = Params.EL0;         % Initial land use emissions
deltaEL     = Params.deltaEL;     % Land use emissions decrease rate

ga          = Params.ga;          % Initial TFP rate
deltaA      = Params.deltaA;      % TFP increase rate

pb          = Params.pb;          % Initial backstop price
deltaPB     = Params.deltaPB;     % Decline rate of backstop price

gsigma      = Params.gsigma;      % Emissions intensity base rate
deltasigma  = Params.deltasigma;  % Decline rate of emissions intensity

f0          = Params.f0;          % Initial forcings of non-CO2 GHGs
f1          = Params.f1;          % Forcings of non-CO2 GHGs in 2100
tforce      = Params.tforce;      % Slope of non-CO2 GHG forcings
        
%% Dynamics / state recursion
i_NEXT = i+1; %time index
F_EX_NEXT = f0 + min(f1-f0, (f1-f0)*(i)/tforce);

% Named functional quantities
theta1 = sigma * (pb/(1000*theta2) * (1-deltaPB)^(i-1));
Gross_Economic_Output = A_TFP*(K^gamma)*((L/1000)^(1-gamma));
Damages = 1 - (a2*(T_AT^a3))/(1+a2*(T_AT^a3)); 
Net_Economic_Output = Damages *(1 - theta1*(mu^theta2))*Gross_Economic_Output;

% ---------------- Climate ------------------------
Radiative_Forcing = xi1*(eta*log((zeta11*M_AT + zeta12*M_UP + xi2*E)/M_AT_Base)/log(2) + F_EX_NEXT);
T_NEXT = Phi_T * T + [Radiative_Forcing; 0];

% ---------------- Carbon Cycle -------------------
M_NEXT = Phi_M * M + [xi2*E; 0; 0];

% ---------------- Economy ------------------------
K_NEXT = (1 - deltaK)^5 * K + 5 * Net_Economic_Output * s;

% ---------------- Auxilliary States --------------
% Time-varying parameters
sigma_NEXT = sigma * exp(-gsigma * (((1-deltasigma)^5)^(i-1)) * 5);
L_NEXT = L * (La/L)^lg;
A_TFP_NEXT = A_TFP / (1 - ga * exp(-deltaA * 5 * (i-1)));
E_LAND_NEXT = EL0*(1-deltaEL)^(i); 

% Intermediate variables
Gross_Economic_Output_NEXT = A_TFP_NEXT*(K_NEXT^gamma)*...
    ((L_NEXT/1000)^(1-gamma));
T_AT_NEXT = T_NEXT(1);
Damages_NEXT = 1 - (a2*(T_AT_NEXT^a3))/(1+a2*(T_AT_NEXT^a3)); 
theta1_NEXT = sigma_NEXT * (pb/(1000*theta2)) * (1-deltaPB)^(i);
Net_Economic_Output_NEXT = Damages_NEXT *(1 - theta1_NEXT*(mu_NEXT^theta2))*...
    Gross_Economic_Output_NEXT;      

% Emission and Consumption (needed for SCC computation)
E_NEXT = 5*(sigma_NEXT*(1-mu_NEXT)*Gross_Economic_Output_NEXT +...
    E_LAND_NEXT);
C_NEXT = 5*(1-s_NEXT) * Net_Economic_Output_NEXT;

% ---------------- Objective ----------------------
J_NEXT = J - L*(((1000/L*C)^(1-alpha) - 1)/((1-alpha)))/(1+rho)^(5*(i-1));          

% numerical scaling
M_NEXT = M_NEXT/sc.M;
L_NEXT = L_NEXT/sc.L;
K_NEXT = K_NEXT/sc.K;
C_NEXT = C_NEXT/sc.C;
J_NEXT = J_NEXT/sc.J;

f = [i_NEXT; T_NEXT; M_NEXT; K_NEXT; ...
     sigma_NEXT; L_NEXT; A_TFP_NEXT; E_LAND_NEXT; F_EX_NEXT; ...
     E_NEXT; C_NEXT; mu_NEXT; s_NEXT; J_NEXT];