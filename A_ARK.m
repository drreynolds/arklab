function Amat = A_ARK(z, Fdata)
% usage: Amat = A_ARK(z, Fdata)
%
% Inputs:  z = current guess for stage solution
%          Fdata = structure containing extra information for evaluating F.
% Outputs: Amat = Jacobian at current guess
%
% This function computes the Jacobian of each intermediate stage residual
% for a multi-stage ARK method, through calling the user-supplied (in
% Fdata) ODE Jacobian function. 
%
% Daniel R. Reynolds
% Department of Mathematics
% Southern Methodist University
% March 2017
% All Rights Reserved

% extract ARK method information from Fdata
Bi = Fdata.Bi;
[Brows, Bcols] = size(Bi);
s   = Bcols - 1;
ci  = Bi(1:s,1);
bi  = (Bi(s+1,2:s+1))';
Ai  = Bi(1:s,2:s+1);
st  = Fdata.stage;
t   = Fdata.t + Fdata.h*ci(st);

% form the ARK Jacobian
Amat = eye(length(z)) - Fdata.h*Ai(st,st)*Fdata.Ji(t, z);

% end of function
