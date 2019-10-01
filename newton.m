function [y,lits,ierr] = newton(Fres, Jres, y0, Fdata, ewt, maxit)
% usage: [y,lits,ierr] = newton(Fres, Jres, y0, Fdata, ewt, maxit)
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
%          ewt = error-weight vector for measuring nonlinear convergence
%          maxit = maximum allowed iterations
% Outputs: y = solution to root-finding problem
%          lits = total # of linear solves taken
%          ierr = output flag denoting success (0) or failure (1)
%
% Note: we use a weighted root-mean squared norm to measure
% solution convergence:
%     ||err|| := sqrt(1/n * sum_{i=1}^n( (err(i)*ewt(i))^2 )) < 1
% where ewt is passed in by the user.
%
% In most cases, the choice 
%     ewt = 1./(rtol*|y0| + atol)
% is appropriate.
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
if (min(ewt) <= 0) 
   error('newton error: illegal error weight vector (all entries must be positive)');
end

% set function to measure convergence (want value <1)
WrmsNorm = @(x) sqrt(sum((x.*ewt).^2)/length(x));

% initialize result, increment vector, residual, statistics
y = y0;
s = ones(size(y));
lits = 0;

% perform iterations
for i=1:maxit

   % compute residual at current guess
   F = Fres(y,Fdata);

   % check residual and increment for stopping
   if (WrmsNorm(s) < 1)
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
