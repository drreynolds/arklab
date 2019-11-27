function udot = fi_timedep_bdry(t, u)
% usage: udot = fi_timedep_bdry(t, u)
%
% Daniel R. Reynolds
% Department of Mathematics
% Southern Methodist University
% October 2019
% All Rights Reserved

% extract problem data
global Pdata;
m  = Pdata.m;
dx = Pdata.dx;
b2 = Pdata.b2;
lambda = Pdata.lambda;

% initialize RHS terms
udot = zeros(m,1);

% interior index set
int = 2:m-1;

% diffusion components
udot(int) = (lambda/dx/dx)*(u(int+1)+u(int-1)-2*u(int));

% end function
