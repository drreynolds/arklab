function u = enforce_timedep_bdry2(t, u)
% usage: u = enforce_timedep_bdry2(t, u)
%
% Daniel R. Reynolds
% Department of Mathematics
% Southern Methodist University
% November 2019
% All Rights Reserved

% extract problem data
global Pdata;
dx = Pdata.dx;
b1 = Pdata.b1;
b2 = Pdata.b2;

% enforce left boundary condition
u(1) = b1(t);

% enforce right boundary condition
u(end) = b2(t);

% end function
