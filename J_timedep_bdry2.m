function J = J_timedep_bdry2(t, u)
% usage: J = J_timedep_bdry2(t, u)
%
% Daniel R. Reynolds
% Department of Mathematics
% Southern Methodist University
% November 2019
% All Rights Reserved

% extract problem data
global Pdata;
m  = Pdata.m;
dx = Pdata.dx;
lambda = Pdata.lambda;

% initialize Jacobian 
J = sparse([],[],[],m,m,3*m);

% diffusion components
dxinv2 = lambda/dx/dx;
for j=2:m-1
   J(j,j-1) = J(j,j-1) + dxinv2;
   J(j,j)   = J(j,j) - 2*dxinv2;
   J(j,j+1) = J(j,j+1) + dxinv2;
end

% end function
