function MPCDICE_test
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
% -------------------------------------------------------------------------
%
%%
clc
close all
display('Testing CasADi installation ...')
input('Press ENTER to continue ...')

try 
    rosenbrock
    display(' ')
    display(' ')
    display(' ')
    display('==============================')
    display('CasADi is working properly.')
    display('You are ready to run MPC-DICE.')
    display('==============================')
catch
    warning('CasADi seems to be not working properly.')
    warning('Please check your installation.')
    try
        import casadi.*
        display('Note: CasADi was found on the MATLAB path.')        
    catch
        warning('Note: CasADi was not found on the MATLAB path.')
        warning('Please check your CasADi installation.')
    end        
end

end

function rosenbrock

import casadi.* 
x1 = SX.sym('x1');
x2 = SX.sym('x2');
x3 = SX.sym('x3');
v = [x1;x2;x3];
f = x1^2 + 100*x3^2;
g = x3 + (1-x1)^2 - x2;
nlp = struct('x', v, 'f', f', 'g', g);

%% Options for IPOPT
opts = struct;
opts.ipopt.max_iter    = 3000;
opts.ipopt.print_level = 5;%0,3
opts.print_time        = 0;
opts.ipopt.acceptable_tol =1e-8;
opts.ipopt.acceptable_obj_change_tol = 1e-8;

solver = nlpsol('solver', 'ipopt', nlp, opts);

res = solver('x0' , [2.5 3.0 0.75],... % solution guess
             'lbx', -inf,...           % lower bound on x
             'ubx',  inf,...           % upper bound on x
             'lbg',    0,...           % lower bound on g
             'ubg',    0);             % upper bound on g
 

end