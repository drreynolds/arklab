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

% right boundary -- enforce via ghost node at m+1
%    b2(t) = u_x(t,pi/2) = (u(m+1)-u(m-1))/(2*dx)
% <=>
%    b2(t)*2*dx + u(m-1) = u(m+1)
% so use "standard" diffusion equation, with this value of the ghost node
ump = b2(t)*2*dx + u(m-1);
udot(m) = (lambda/dx/dx)*(ump+u(m-1)-2*u(m));

% end function
