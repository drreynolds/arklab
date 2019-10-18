function J = J_timedep_bdry(t, u)
% usage: J = J_timedep_bdry(t, u)
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

% initialize Jacobian 
J = sparse([],[],[],m,m,3*m);

% diffusion components
dxinv2 = 1/dx/dx;
for j=2:m-1
   J(j,j-1) = J(j,j-1) + dxinv2;
   J(j,j)   = J(j,j) - 2*dxinv2;
   J(j,j+1) = J(j,j+1) + dxinv2;
end
J(m,m-1) = 2*dxinv2;
J(m,m) = -2*dxinv2;

% end function
