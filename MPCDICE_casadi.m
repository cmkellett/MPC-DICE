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
% This code provides a receding horizon (MPC) implementation of DICE. 
% The model implemented here differs slightly from that of the GAMS code
% available from Prof. W. Nordhaus at
% (http://aida.wss.yale.edu/~nordhaus/homepage/Web-DICE-2013-April.htm).
% Using the default parameters, the numerical results, particularly for the
% social cost of carbon, are very close.  However, more extreme parameter
% values are known to give rise to problems in the GAMS implementation.
% Such issues should not arise with this code.  See the accompanying 
% documentation for more details.
%
% -------------------------------------------------------------------------
% Release notes
%
%  - 1.0 Initial release
% -------------------------------------------------------------------------

clc
clear all
close all

% import casadi name space
import casadi.*

%% ========================================================================
% Define data of the DICE MPC-loop & construct NLP
% =========================================================================

% Set parameters using an appropriate function
Params = assign_parameters_v2016;

Data.N = Params.N;  % Prediction horizon set in above script file
Data.step = 1;      % 1 == 1-step MPC strategy
Data.t0 = 1;        % initial simulation time; needs to be >=1
Data.tf = 10;       % final simulation time; set tf = 1 to solve single OCP
Data.nx = 1+6+9+1;  % total # of states: time; 6 "endogenous" states; 
                    %   9 auxilliary states including 5 "exogenous" states,
                    %   consumption, emissions, and shifted inputs;
                    %   objective
Data.nu = 2;        % # of inputs


Params.x0 = set_initial_conditions(Data.t0,Params); % get x0

display('============================================================')
display('                        MPC-DICE')
display('                   Copyright (C) 2018 by')
display('Timm Faulwasser, Christopher M. Kellett and Steven R. Weller')
display(' ')
display(' ')
display('============================================================')
%% Symbolic problem construction
display(['Symbolic construction of DICE optimal control problem ...'])
tic
[SingleOCP, SequenceOCP] = ConstructNLP(Data.N, Params.x0, Params); %construct NLPs
solver_ini = nlpsol('solver', 'ipopt', SingleOCP.nlp, Params.opts); % create IPOPT solver object    
solver     = nlpsol('solver', 'ipopt', SequenceOCP.nlp, Params.opts); % create IPOPT solver object
t_NLP = toc;

% Prepare output data struct
Data.state_list = {'time', 'T_AT', 'T_LO', 'M_AT', 'M_UP',...
                   'M_LO', 'K', 'sigma', 'L', 'A_TFP', 'Emission', 'C', ...
                   'mu', 's', 'J'};               
Data.xMPC  = [];
Data.uMPC  = [];
Data.wMPC  = [];
Data.lam_C = [];
Data.lam_E = [];
Data.t_SOL = [];

%% Main loop
switch Data.tf
    case 1 % Solve single OCP
                             
        %solve the NLP 
        tic       
        res = solver_ini( 'x0' , SingleOCP.xu_guess,...  % solution guess (hot start)
                          'lbx', SingleOCP.xu_LB,...     % lower bound on x
                          'ubx', SingleOCP.xu_UB,...     % upper bound on x
                          'lbg', 0, ...     % lower bound on g
                          'ubg', 0);        % parameters of NLP 
        Data.t_SOL = [Data.t_SOL; toc];

        xu_opt = full(res.x);
        index_w = (Data.N+1)*Data.nu -1;
        index_x = (Data.N+1)*Data.nx;
        w_opt = xu_opt(end-index_w:end);        
        w_opt = reshape(w_opt, Data.nu, (Data.N+1));
        
        x_opt = xu_opt(1:index_x);
        x_opt = reshape(x_opt, Data.nx, (Data.N+1));
        u_opt = x_opt(end-2:end-1, 1:end);
        
        Data.xMPC = x_opt(:, 1:Data.N+1);
        Data.uMPC = u_opt;
        Data.wMPC = w_opt;
        
        lambda  = full(res.lam_g);
        lambda = lambda(1:end-15);
        
        lambda = reshape(lambda, Data.nx, Data.N+1);
        Data.lam_E = lambda(13,:); % dW/dE
        Data.lam_C = lambda(14,:)./Params.sc.C; % dW/dC
        Data.SCC   = -1000*Data.lam_E./Data.lam_C;
        
        % print computation times to screen
        display('============================================================')
        display(['Time to construct NLP:    ', num2str(t_NLP)])
        display(['Time to solve single NLP: ', num2str(max(Data.t_SOL))])
        display('============================================================')
            
    otherwise %% Run MPC loop           
        
        for k = Data.t0 : Data.step : Data.tf
            display(['MPC step k = ', num2str(k), ' of ' num2str(ceil(Data.tf/Data.step))])

            % assign initial condition
            if k > Data.t0
                xk = Data.xMPC(:,end); 
            end

            % solve the NLP
            tic
            if k == Data.t0 % very first OCP
                res = solver_ini( 'x0' , SingleOCP.xu_guess,...% solution guess (hot start)
                                  'lbx', SingleOCP.xu_LB,...   % lower bound on x
                                  'ubx', SingleOCP.xu_UB,...   % upper bound on x
                                  'lbg', 0,...                 % lower bound on g
                                  'ubg', 0);                   % parameters of NLP                 
            else  % sequence OCPs                  
                res = solver( 'x0' , xu_guess,...              % solution guess (hot start)
                              'lbx', SequenceOCP.xu_LB,...     % lower bound on x
                              'ubx', SequenceOCP.xu_UB,...     % upper bound on x
                              'lbg', 0,...                     % lower bound on g
                              'ubg', 0,...                     % upper bound on g
                              'p',   xk);                      % parameters of NLP                                                
            end 
            
            xu_opt = full(res.x);
            index_w = (Data.N+1)*Data.nu -1;
            index_x = (Data.N+1)*Data.nx;
            w_opt = xu_opt(end-index_w:end);        
            w_opt = reshape(w_opt, Data.nu, (Data.N+1)); 
            x_opt = xu_opt(1:index_x);
            x_opt = reshape(x_opt, Data.nx, (Data.N+1));
            
            if k == Data.t0
                Data.xMPC = x_opt(:,1);
            end
                        
            u_opt = x_opt(end-2:end-1, [1:1:Data.step]); %actual input
            x_opt = x_opt(:, 1+[1:1:Data.step]); %state
            w_opt = w_opt(:, [1:1:Data.step]); % shifted input = optimization input
            
            tic
            Data.t_SOL = [Data.t_SOL; toc]; 
            Data.xMPC = [Data.xMPC, x_opt];
            Data.uMPC = [Data.uMPC, u_opt]; 
            Data.wMPC = [Data.wMPC, w_opt]; 
            
            lambda  = full(res.lam_g);
            if k == Data.t0
                lambda = lambda(1:end-15);
                lambda = reshape(lambda, Data.nx, Data.N+1);
            else
                lambda = reshape(lambda, Data.nx, Data.N+1);
            end
            
            lambdaE = lambda(13,1:1:Data.step); %11
            lambdaC = lambda(14,1:1:Data.step)/Params.sc.C; %12
            
            Data.lam_C = [Data.lam_C, lambdaC];
            Data.lam_E = [Data.lam_E, lambdaE];

            % construct guess
            if k <= Data.tf
                    if k == 1                    
                        offset = (Data.N+1)*Data.nu + 1*Data.nx;
                        xu_guess = [xu_opt(Data.nx+length(x_opt(:))+1:end-offset); x_opt(:); xu_opt(end-offset+length(u_opt(:))+1:end); u_opt(:)];
                    else
                        offset = (Data.N+1)*Data.nu;
                        xu_guess = [xu_opt(length(x_opt(:))+1:end-offset); x_opt(:); xu_opt(end-offset+length(u_opt(:))+1:end); u_opt(:)];
                    end
            end

        end   
    
    % print computation times to screen
    display('============================================================')
    display(['Time to construct NLP:    ', num2str(t_NLP)])
    display(['Minimal time to solve NLP: ', num2str(min(Data.t_SOL))])
    display(['Maximal time to solve NLP: ', num2str(max(Data.t_SOL))])
    display(['Average time to solve NLP: ', num2str(mean(Data.t_SOL))])
    display(['Number of NLPs solved:     ', num2str(ceil(Data.tf/Data.step))])
    display('============================================================')
end

%% User-friendly unpacking of results

%% DICE "endogenous" states
TATM = Data.xMPC(2,:);
TLO = Data.xMPC(3,:);
MATM = Params.sc.M*Data.xMPC(4,:);
MUP = Params.sc.M*Data.xMPC(5,:);
MLO = Params.sc.M*Data.xMPC(6,:);
K = Params.sc.K*Data.xMPC(7,:);

%% DICE "exogenous" states
sigma = Data.xMPC(8,:);
L = Params.sc.L*Data.xMPC(9,:);
A_TFP = Data.xMPC(10,:);

Emissions = Data.xMPC(13,:);  
Consumption = Params.sc.C*Data.xMPC(14,:);

%% DICE inputs
Savings_Rate = Data.uMPC(2,:);
miu = Data.uMPC(1,:);

%% SC-CO2
Data.SCC   = -1000*Data.lam_E./Data.lam_C;
SCC = Data.SCC;

if Data.tf == 1
    end_year = Params.BaseYear+5*Data.N;
    horizon = Data.N+1;
else
    end_year = Params.BaseYear+5*Data.tf;
    horizon = Data.tf;
end
years = Params.BaseYear:5:end_year;

%% Compute auxiliary outputs
for i=1:horizon
    Gross_Economic_Output(i) = A_TFP(i)*(K(i)^Params.gamma)*...
        ((L(i)/1000)^(1-Params.gamma));
    Damages_Fraction(i) = (Params.a2*(TATM(i)^Params.a3))/...
        (1+Params.a2*(TATM(i)^Params.a3)); 
    Abatement_Fraction(i) = (Params.pb/(1000*Params.theta2))*...
        ((1-Params.deltaPB)^(i-1))*sigma(i)*(miu(i)^Params.theta2);
    Net_Output(i) = (1-Damages_Fraction(i))*(1 - Abatement_Fraction(i))*...
        Gross_Economic_Output(i);    
    E_Industrial(i) = sigma(i)*(1-miu(i))*Gross_Economic_Output(i);
    Per_Cap_Consumption(i) = 1000*(1-Savings_Rate(i)) * Net_Output(i) / L(i);
    Atm_Carbon_ppm(i) = MATM(i)/2.13;
    Marg_Cost_Abatement(i) = Params.pb*((1-Params.deltaPB)^(i-1))*... 
        (miu(i)^(Params.theta2-1));
end



%% Save results in unique file
SaveData.Data = Data;
SaveData.BaseYear = Params.BaseYear;
SaveData.end_year = end_year;
SaveData.years = years;
SaveData.TATM = TATM;
SaveData.TLO = TLO;
SaveData.MATM = MATM;
SaveData.MUP = MUP;
SaveData.MLO = MLO;
SaveData.K = K;
SaveData.sigma = sigma;
SaveData.L = L;
SaveData.A_TFP = A_TFP;
SaveData.Emissions = Emissions;
SaveData.Consumption = Consumption;
SaveData.miu = miu;
SaveData.Savings_Rate = Savings_Rate;
SaveData.SCC = SCC;
SaveData.Per_Cap_Consumption = Per_Cap_Consumption;
SaveData.Gross_Economic_Output = Gross_Economic_Output;

file_id = sprintf('%s%s%s','DataResults_',Params.parameter_set,'.mat');
save(file_id,'SaveData');

%% Sample Plotting

figure(1)
stem(years(1:length(years)-1),SCC), axis([Params.BaseYear end_year 0 max(SCC)+100])
title('Social Cost of Carbon')
xlabel('Years')

figure(2)
subplot(2,1,1), stem(years(1:length(years)-1),miu), axis([Params.BaseYear end_year 0 1.2])
ylabel('miu')
subplot(2,1,2), stem(years(1:length(years)-1),Savings_Rate), axis([Params.BaseYear end_year 0 0.5])
ylabel('Savings Rate')

figure(3)
stem(years(1:length(years)-1),Per_Cap_Consumption), axis([Params.BaseYear end_year 0 max(Per_Cap_Consumption)+5])
title('Per Capita Consumption')
xlabel('Years')

figure(4)
stem(years,TATM)
title('Atmospheric Temperature Anomaly')
xlabel('Years')

% Tidy Workspace.  Comment out the following two commands to maintain all
%   variables in workspace.
clear_list = {'index_w','index_x','lambda','lambda_G','res','SequenceOCP','SingleOCP',...
    'solver','solver_ini','t_NLP','u_opt','w_opt','x_opt','xu_opt'};
clear(clear_list{:});
clear('clear_list');