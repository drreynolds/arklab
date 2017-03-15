function dy = fi_p3(t, y)
% usage: dy = fi_p3(t, y)
%
% Daniel R. Reynolds
% Department of Mathematics
% Southern Methodist University
% March 2017
% All Rights Reserved

% extract problem data
global Pdata;
d1 = Pdata.d1;
d2 = Pdata.d2;
m  = Pdata.m;
dx = Pdata.dx;

% extract solution components
u = y(1:m);
v = y(m+1:2*m);

% initialize RHS terms
du = zeros(m,1);
dv = zeros(m,1);

% enforce stationary boundary conditions
du(1) = 0;  du(m) = 0;  dv(1) = 0;  dv(m) = 0;

% diffusion components
du(2:m-1) = d1/dx/dx*(u(3:m)+u(1:m-2)-2*u(2:m-1));
dv(2:m-1) = d2/dx/dx*(v(3:m)+v(1:m-2)-2*v(2:m-1));

% combine together into output
dy = [du; dv];

% end function
