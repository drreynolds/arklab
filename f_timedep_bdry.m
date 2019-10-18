function udot = f_timedep_bdry(t, u)
% usage: udot = f_timedep_bdry(t, u)
%
% Daniel R. Reynolds
% Department of Mathematics
% Southern Methodist University
% October 2019
% All Rights Reserved

% extract problem data
global Pdata;
c  = Pdata.c;
m  = Pdata.m;
dx = Pdata.dx;
b2 = Pdata.b2;
xspan = Pdata.xspan;
b1dot = Pdata.b1dot;

% initialize RHS terms
udot = zeros(m,1);

% interior index set
int = 2:m-1;
x = xspan(int);

% diffusion components
udot(int) = (1/dx/dx)*(u(int+1)+u(int-1)-2*u(int));

% forcing term
udot(int) = udot(int) + c(1)*cos(x*1) + c(2)*cos(x*2) + c(3)*cos(x*3) + c(4)*cos(x*4) + c(5)*cos(x*5);

% left boundary
udot(1) = b1dot(t);

% right boundary -- enforce via ghost node at m+1
%    b2(t) = u_x(t,pi/2) = (u(m+1)-u(m-1))/(2*dx)
% <=>
%    b2(t)*2*dx + u(m-1) = u(m+1)
% so use "standard" equation, with this value of the ghost node
ump = b2(t)*2*dx + u(m-1);
udot(m) = (1/dx/dx)*(ump+u(m-1)-2*u(m)) - c(2) + c(4);

% end function
