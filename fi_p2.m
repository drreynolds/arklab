function dy = fi_p2(t, y)
% usage: dy = fi_p2(t, y)
%
% Daniel R. Reynolds
% Department of Mathematics
% Southern Methodist University
% March 2017
% All Rights Reserved

% model parameters
global Pdata;
b  = Pdata.b; 
ep = Pdata.ep;

% form the ODE RHS
dy = [0;
      0;
      (b-y(3))/ep];

% end function
