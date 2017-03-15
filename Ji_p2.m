function J = Ji_p2(t, y)
% usage: J = Ji_p2(t, y)
%
% Daniel R. Reynolds
% Department of Mathematics
% Southern Methodist University
% March 2017
% All Rights Reserved

% model parameters
global Pdata;
ep = Pdata.ep;

% form the ODE Jacobian
J = [0, 0, 0;
     0, 0, 0;
     0, 0, -1/ep];

% end function
