function dy = fe_p1(t, y)
% usage: dy = fe_p1(t, y)
%
% Daniel R. Reynolds
% Department of Mathematics
% Southern Methodist University
% March 2017
% All Rights Reserved

% extract variables
u = y(1);
v = y(2);

% form the ODE RHS
du = v;
dv = -u;
dy = [du; dv];

% end function
