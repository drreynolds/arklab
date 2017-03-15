function dy = fe_p2(t, y)
% usage: dy = fe_p2(t, y)
%
% Daniel R. Reynolds
% Department of Mathematics
% Southern Methodist University
% March 2017
% All Rights Reserved

% model parameters
global Pdata;
a  = Pdata.a; 

% form the ODE RHS
dy = [a - (y(3)+1)*y(1) + y(1)*y(1)*y(2);
      y(3)*y(1) - y(1)*y(1)*y(2);
      -y(3)*y(1)];

% end function
