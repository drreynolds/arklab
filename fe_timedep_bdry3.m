function udot = fe_timedep_bdry3(t, u)
% usage: udot = fe_timedep_bdry3(t, u)
%
% Daniel R. Reynolds
% Department of Mathematics
% Southern Methodist University
% October 2019
% All Rights Reserved

% extract problem data
global Pdata;
c = Pdata.c;
m = Pdata.m;
xspan = Pdata.xspan;

% initialize RHS terms
udot = zeros(m,1);

% interior index set
int = 2:m-1;
x = xspan(int);

% forcing term
udot(int) = c(1)*cos(x*1) + c(2)*cos(x*2) + c(3)*cos(x*3) + c(4)*cos(x*4) + c(5)*cos(x*5);

% handle time-dependent boundary condition approximations
global bDotApprox

% on the first call, evaluate time-dependent boundary values to fill approximation structure
if (bDotApprox.nstored == 0)
  delta = bDotApprox.h/bDotApprox.nmax;
  tvals = linspace(t-bDotApprox.nmax/2*delta, t+bDotApprox.nmax/2*delta, bDotApprox.nmax);
  for i=1:bDotApprox.nmax
    bDotApprox.b1(i) = Pdata.b1(tvals(i));
    bDotApprox.b2(i) = Pdata.b2(tvals(i));
    bDotApprox.t(i) = tvals(i);
  end
  bDotApprox.nstored = bDotApprox.nmax;

% on other calls, only update approximation structure at new "t" value
else
  if (min(abs(t-bDotApprox.t)) > 10*sqrt(eps))
    [~,i] = max(abs(t-bDotApprox.t));
    bDotApprox.b1(i) = Pdata.b1(t);
    bDotApprox.b2(i) = Pdata.b2(t);
    bDotApprox.t(i) = t;
  end
end

% left boundary
udot(1) = interp_deriv(bDotApprox.t, bDotApprox.b1, t);

% right boundary
udot(m) = interp_deriv(bDotApprox.t, bDotApprox.b2, t);



function dp = interp_deriv(x,y,t)
% This routine evaluates the first-derivative of the
% Newton interpolating polynomial, p'(t).
%
% Inputs:   x   nodal locations [array of length n] 
%               (assume x(i) ~= x(j) for i ~= j)
%           y   data values [array of length n]
%           t   evaluation point
% Output:   dp  p'(t)

% compute Newton coefficients using divided differences; 
% store in matrix d (diagonal holds Newton coefficients)
n = length(x);
d = zeros(n, n);              % initialize d
d(:,1) = reshape(y, n, 1);    % fill first column with contents of y
for j=2:n
   for i=j:n                  % perform divided-differences algorithm
      d(i,j) = (d(i-1,j-1) - d(i,j-1)) / (x(i-j+1) - x(i));
   end
end

% evaluate reusable factors
dt = t-x;

% evaluate Newton polynomial derivative
dp = 0;
for i=1:n
  dphi = 0;
  for k=1:i-1
    tprod = 1;
    for j=1:i-1
      if (j == k), continue, end
      tprod = tprod * dt(j);
    end
    dphi = dphi + tprod;
  end
  dp = dp + d(i,i)*dphi;
end

% end function
