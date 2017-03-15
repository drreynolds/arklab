function dy = fe_p3(t, y)
% usage: dy = fe_p3(t, y)
%
% Daniel R. Reynolds
% Department of Mathematics
% Southern Methodist University
% March 2017
% All Rights Reserved

% extract problem data
global Pdata;
a  = Pdata.a; 
b  = Pdata.b; 
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

% reaction components
du(2:m-1) = a - (b+1)*u(2:m-1) + u(2:m-1).*u(2:m-1).*v(2:m-1);
dv(2:m-1) = b*u(2:m-1) - u(2:m-1).*u(2:m-1).*v(2:m-1);

% combine together into output
dy = [du; dv];

% end function
