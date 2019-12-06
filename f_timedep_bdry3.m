function udot = f_timedep_bdry3(t, u)
% usage: udot = f_timedep_bdry3(t, u)
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
lambda = Pdata.lambda;
xspan = Pdata.xspan;

% initialize RHS terms
udot = zeros(m,1);

% interior index set
int = 2:m-1;
x = xspan(int);

% diffusion components
udot(int) = (lambda/dx/dx)*(u(int+1)+u(int-1)-2*u(int));

% forcing term
udot(int) = udot(int) + c(1)*cos(x*1) + c(2)*cos(x*2) + c(3)*cos(x*3) + c(4)*cos(x*4) + c(5)*cos(x*5);


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
p = polyfit(bDotApprox.t, bDotApprox.b1, bDotApprox.nmax-1);
dp = polyder(p);
udot(1) = polyval(dp,t);

% right boundary
p = polyfit(bDotApprox.t, bDotApprox.b2, bDotApprox.nmax-1);
dp = polyder(p);
udot(m) = polyval(dp,t);

% end function
