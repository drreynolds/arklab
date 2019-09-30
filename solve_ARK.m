function [tvals, Y, nsteps, lits] = solve_ARK(fe, fi, Ji, tvals, Y0, Be, ...
                                              Bi, rtol, atol, hmin, hmax)
% usage: [tvals, Y, nsteps, lits] = solve_ARK(fe, fi, Ji, tvals, Y0, Be, ...
%                                             Bi, rtol, atol, hmin, hmax)
%
% Adaptive time step additive Runge-Kutta solver for the
% vector-valued ODE problem  
%     y' = fe(t,Y) + fi(t,Y), t in tvals, y in R^m,
%     Y(t0) = [y1(t0), y2(t0), ..., ym(t0)]'.
%
% Inputs:
%     fe     = function handle for fe(t,Y)
%     fi     = function handle for fi(t,Y)
%     Ji     = function handle for Jacobian of fi, J(t,Y)
%     tvals  = [t0, t1, t2, ..., tN]
%     Y0     = initial value array (column vector of length m)
%     Be,Bi  = Butcher table matrices for ARK coefficients, of the form
%                 Be = [ce Ae;       Bi = [ci Ai;
%                       q  be;             q  bi;
%                       p  be2 ]           p  bi2 ]
%              Here, ce,ci are vectors of stage time fractions (s-by-1),
%                    Ae,Ai are matrices of Butcher coefficients (s-by-s),
%                    q is an integer denoting the method order of accuracy,
%                    be,bi are vectors of solution weights (1-by-s),
%                    p is an integer denoting the embedding order of accuracy,
%                    be2,bi2 are vectors of embedding weights (1-by-s),
%              The [p, be2] and [p, bi2] rows are optional.  If
%              both of those are not provided the method will default to
%              taking fixed step sizes of size hmin.
%     rtol   = desired relative error of solution  (scalar)
%     atol   = desired absolute error of solution  (vector or scalar)
%     hmin   = minimum internal time step size (hmin <= t(i)-t(i-1), for all i)
%     hmax   = maximum internal time step size (hmax >= hmin)
%
% Outputs: 
%     tvals  = the same as the input array tvals
%     y      = [y(t0), y(t1), y(t2), ..., y(tN)], where each
%               y(t*) is a column vector of length m.
%     nsteps = number of internal time steps taken by method
%     lits   = number of linear solves required by method
%
% Note: to run in fixed-step mode, call with hmin=hmax as the desired 
% time step size, and set the tolerances to large positive numbers.
%
% Daniel R. Reynolds
% Department of Mathematics
% Southern Methodist University
% March 2017
% All Rights Reserved

% check for compatible Be,Bi tables
if (size(Be,2) ~= size(Bi,2))   
   error('solve_ARK error: Be and Bi must have the same number of stages')
end
s = size(Be,2) - 1;          % number of stages
if (Be(s+1,1) ~= Bi(s+1,1))
   error('solve_ARK error: Be and Bi must have the same method order')
end
if (size(Be,1) > size(Be,2))
   if (Be(s+1,2) ~= Bi(s+1,2))
      error('solve_ARK error: Be and Bi must have the same embedding order')
   end
end

% extract ERK method information from Be
[Brows, Bcols] = size(Be);
ce = Be(1:s,1);         % stage time fraction array
be = (Be(s+1,2:s+1))';  % solution weights (convert to column)
Ae = Be(1:s,2:s+1);     % RK coefficients
q  = Be(s+1,1);         % method order

% extract DIRK method information from Bi
[Brows, Bcols] = size(Bi);
ci = Bi(1:s,1);         % stage time fraction array
bi = (Bi(s+1,2:s+1))';  % solution weights (convert to column)
Ai = Bi(1:s,2:s+1);     % RK coefficients

% initialize as non-embedded, until proven otherwise
embedded = 0;
p = 0;
if (Brows > Bcols)
   if ((max(abs(Be(s+2,2:s+1))) > eps) && (max(abs(Bi(s+2,2:s+1))) > eps))
      embedded = 1;
      be2 = (Be(s+2,2:s+1))';
      bi2 = (Bi(s+2,2:s+1))';
      p = Be(s+2,1);
   end
end

% initialize output arrays
N = length(tvals);
m = length(Y0);
Y = zeros(m,N);
Y(:,1) = Y0;

% initialize diagnostics
c_fails = 0;   % total convergence failures
a_fails = 0;   % total accuracy failures

% set the solver parameters
newt_maxit = 20;           % max number of Newton iterations
newt_ftol  = 1e-10;        % Newton solver residual tolerance
newt_stol  = 1e-10;        % Newton solver solution tolerance
h_reduce   = 0.1;          % failed step reduction factor 
h_safety   = 0.9;          % adaptivity safety factor
h_growth   = 10;           % adaptivity growth bound
ONEMSM     = 1-sqrt(eps);  % coefficients to account for
ONEPSM     = 1+sqrt(eps);  %   floating-point roundoff
ERRTOL     = 1.1;          % upper bound on allowed step error
                           %   (in WRMS norm)

% initialize temporary variables
t = tvals(1);
Ynew = Y0;

% create Fdata structure for Newton solver and step solutions
Fdata.fe = fe;    % ODE RHS function names
Fdata.fi = fi;
Fdata.Ji = Ji;    % ODE RHS Jacobian function name
Fdata.Be = Be;    % Butcher tables
Fdata.Bi = Bi;
Fdata.s  = s;     % number of stages

% set function names for Newton solver residual/Jacobian
Fun = @F_ARK;
Jac = @A_ARK;

% set initial time step size
h = hmin;

% initialize work counters
nsteps = 0;
lits   = 0;

% iterate over output time steps
for tstep = 2:length(tvals)

   % loop over internal time steps to get to desired output time
   while (t < tvals(tstep)*ONEMSM)
      
      % bound internal time step 
      h = max([h, hmin]);            % enforce minimum time step size
      h = min([h, hmax]);            % maximum time step size
      h = min([h, tvals(tstep)-t]);  % stop at output time

      % set Fdata values for this step
      Fdata.h    = h;    % current step size
      Fdata.yold = Y0;   % solution from previous step
      Fdata.t    = t;    % time of last successful step

      % initialize data storage for multiple stages
      z = zeros(m,s);

      % reset stage failure flag
      st_fail = 0;
      
      % loop over stages
      for stage = 1:s
         
         % set Newton initial guess as previous stage solution
         Yguess = Ynew;
         
         % set current stage index into Fdata structure
         Fdata.stage = stage;
         
         % construct RHS comprised of old time data
         %    zi = y_n + h*sum_{j=1}^s (a(i,j)*fj)
         % <=>
         %    zi - h*(a(i,i)*fi) = y_n + h*sum_{j=1}^{i-1} (a(i,j)*fj)
         % =>
         %    rhs = y_n + h*sum_{j=1}^{i-1} (a(i,j)*fj)
         Fdata.rhs = Y0;
         for j = 1:stage-1
            Fdata.rhs = Fdata.rhs + h*Ae(stage,j)*fe(t+h*ce(j), z(:,j)) ...
                                  + h*Ai(stage,j)*fi(t+h*ci(j), z(:,j));
         end
         
         % call Newton solver to compute new stage solution
         [Ynew,lin,ierr] = newton(Fun, Jac, Yguess, Fdata, ...
                                  newt_ftol, newt_stol, newt_maxit);

         % increment total linear solver statistics
         lits = lits + lin;
         
         % if Newton method failed, set relevant flags/statistics
         % and break out of stage loop
         if (ierr ~= 0) 
            st_fail = 1;
            c_fails = c_fails + 1;
            break;
         end
         
         % store stage solution
         z(:,stage) = Ynew;
         
      end
      
      % increment number of internal time steps taken
      nsteps = nsteps + 1;
      
      % compute new solution (and embedding if available)
      [Ynew,Y2] = Y_ARK(z,Fdata);

      % if stages succeeded and time step adaptivity enabled, check step accuracy
      if ((st_fail == 0) & embedded)

         % estimate error in current step
         err_step = max(norm((Ynew - Y2)./(rtol*Ynew + atol),inf), eps);
         
         % if error too high, flag step as a failure (will be be recomputed)
         if (err_step > ERRTOL*ONEPSM) 
            a_fails = a_fails + 1;
            st_fail = 1;
         end
         
      end

      % if step was successful (solves succeeded, and error acceptable)
      if (st_fail == 0) 
         
         % update solution and time for last successful step
         Y0 = Ynew;
         t  = t + h;
         
         % for embedded methods, use error estimate to adapt the time step
         if (embedded) 

            h_old = h;
            if (err_step == 0.0)     % no error, set max possible
               h = tvals(end)-t;
            else                     % set next h (I-controller)
               h = h_safety * h_old * err_step^(-1.0/q);
            end

            % enforce maximum growth rate on step sizes
            h = min(h_growth*h_old, h);

         % otherwise, just use the fixed minimum input step size
         else
            h = hmin;
         end
         
      % if step solves or error test failed
      else

         % if already at minimum step, just return with failure
         if (h <= hmin) 
            fprintf('Cannot achieve desired accuracy.\n');
            fprintf('Consider reducing hmin or increasing rtol.\n');
            return
         end

         % otherwise, reset guess, reduce time step, retry solve
         Ynew = Y0;
         h = h * h_reduce;
         
      end  % end logic tests for step success/failure
      
   end  % end while loop attempting to solve steps to next output time

   % store updated solution in output array
   Y(:,tstep) = Ynew;
   
end  % time step loop

% end solve_ARK function
end




function [y,y2] = Y_ARK(z, Fdata)
% usage: [y,y2] = Y_ARK(z, Fdata)
%
% Inputs:
%    z     = stage solutions [z1, ..., zs]
%    Fdata = structure containing extra problem information
%
% Outputs: 
%    y     = step solution built from the z values
%    y2    = embedded solution (if embedding included in Butcher 
%               table; otherwise the same as y)

% extract method information from Fdata
Be = Fdata.Be;
Bi = Fdata.Bi;
[Brows, Bcols] = size(Be);
s = Bcols - 1;
ce = Be(1:s,1);
ci = Bi(1:s,1);
be = (Be(s+1,2:s+1))';
bi = (Bi(s+1,2:s+1))';

% check to see if we have coefficients for embedding
if (Brows > Bcols)
   be2 = (Be(s+2,2:s+1))';
   bi2 = (Bi(s+2,2:s+1))';
else
   be2 = be;
   bi2 = bi;
end

% get some problem information
[zrows,zcols] = size(z);
nvar = zrows;
if (zcols ~= s)
   error('Y_ARK error: z has incorrect number of stages');
end

% call RHS at our stages
fe = zeros(nvar,s);
fi = zeros(nvar,s);
for is=1:s
   t = Fdata.t + Fdata.h*ce(is);
   fe(:,is) = Fdata.fe(t, z(:,is));

   t = Fdata.t + Fdata.h*ci(is);
   fi(:,is) = Fdata.fi(t, z(:,is));
end

% form the solutions
%    ynew = yold + h*sum(be(j)*fe(j) + bi(j)*fi(j))
y  = Fdata.yold + Fdata.h*fe*be + Fdata.h*fi*bi;
y2 = Fdata.yold + Fdata.h*fe*be2 + Fdata.h*fi*bi2;

% end of function
end




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
end




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
end
