function u = enforce_timedep_bdry(t, u)
% usage: u = enforce_timedep_bdry(t, u)
%
% Daniel R. Reynolds
% Department of Mathematics
% Southern Methodist University
% October 2019
% All Rights Reserved

% extract problem data
global Pdata;
dx = Pdata.dx;
b1 = Pdata.b1;
b2 = Pdata.b2;

% enforce left boundary condition
u(1) = b1(t);

% end function
