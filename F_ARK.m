function F = F_ARK(z, Fdata)
% usage: F = F_ARK(z, Fdata)
%
% Inputs:  z = current guess for stage solution
%          Fdata = structure containing extra information for evaluating F.
% Outputs: F = residual at current guess
%
% This function computes the (non)linear residuals for an intermediate
% stage solution, through calling the user-supplied (in Fdata) ODE
% right-hand side function.
%
% Daniel R. Reynolds
% Department of Mathematics
% Southern Methodist University
% March 2017
% All Rights Reserved

% extract ARK method information from Fdata
Bi = Fdata.Bi;
[Brows, Bcols] = size(Bi);
s  = Bcols - 1;
ci = Bi(1:s,1);
Ai = Bi(1:s,2:s+1);
h  = Fdata.h;
st = Fdata.stage;
t  = Fdata.t + Fdata.h*ci(st);

% form the ARK residual
%    F = z - rhs - h*(ai(stage,stage)*fstage)
F = z - Fdata.rhs - h*Ai(st,st)*Fdata.fi(t, z);

% end of function
