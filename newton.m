function [y,lits,ierr] = newton(Fres, Jres, y0, Fdata, rtol, atol, maxit)
% usage: [y,lits,ierr] = newton(Fres, Jres, y0, Fdata, rtol, atol, maxit)
%
% Newton solver for the root-finding problem defined by the function Fres,
%     Fres(y,Fdata) = 0
%
% Inputs:  Fres = function name for nonlinear residual, F(y,Fdata).  Note, all 
%                data required to evaluate F (other than y) should be stored
%                in the data structure Fdata.
%          Jres = function name for Jacobian of nonlinear residual, 
%                A = partial_y F(y,Fdata).  
%                Jres should use the same data structure for additional data
%                as F.
%          y0 = initial guess
%          Fdata = structure containing extra information for evaluating F.
%          rtol = desired relative solution tolerance
%          atol = desired absolute solution tolerance
%          maxit = maximum allowed iterations
% Outputs: y = solution to root-finding problem
%          lits = total # of linear solves taken
%          ierr = output flag denoting success (0) or failure (1)
%
% Note: we use the inputs rtol and atol to define a weighted
% root-mean squared norm to measure solution convergence:
%     ||err|| := sqrt(1/n * sum_{i=1}^n( (err(i)*w(i))^2 ))
% where  w = 1./(rtol*|y0| + atol).
%
% Daniel R. Reynolds
% Department of Mathematics
% Southern Methodist University
% October 2019
% All Rights Reserved

% check solver inputs
if (maxit < 1) 
   error('newton error: requires at least 1 iteration (maxit)');
end
if (rtol < 0) 
   error('newton error: relative tolerance must be non-negative (rtol)');
end
if (atol < 0) 
   error('newton error: absolute tolerance must be non-negative (atol)');
end
if ((rtol == 0) && (atol == 0))
   error('newton error: at least one tolerance must be positive (rtol,atol)');
end

% set function to measure convergence (want value <1)
w = 1./(rtol*abs(y0)+atol);
wrmsnorm = @(x) sqrt(sum((x.*w).^2)/length(x));

% initialize result, increment vector, residual, statistics
y = y0;
s = ones(size(y));
lits = 0;

% perform iterations
for i=1:maxit

   % compute residual at current guess
   F = Fres(y,Fdata);

   % check residual and increment for stopping
   if (wrmsnorm(s) < 1)
      ierr = 0;
      return
   end
   
   % compute Jacobian
   A = Jres(y,Fdata);
   
   % perform Newton update
   s = A\F;
   y = y - s;
   lits = lits + 1;
   
end

% if we've made it to this point, the Newton iteration did not converge
ierr = 1;
fprintf('\nnewton warning: nonconvergence after %i iterations (|F| = %g)\n',maxit,norm(F,inf));

% end of function
