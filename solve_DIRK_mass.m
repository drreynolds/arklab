function [tvals,Y,nsteps,lits] = solve_DIRK_mass(Mn,fcn,Jfcn,tvals,Y0,B,rtol,atol,hmin,hmax,alg)
% usage: [tvals,Y,nsteps,lits] = solve_DIRK_mass(Mn,fcn,Jfcn,tvals,Y0,B,rtol,atol,hmin,hmax,alg)
%
% Adaptive time step diagonally-implicit Runge-Kutta solver for the
% vector-valued ODE problem 
%     y' = F(t,Y), t in tvals, y in R^m,
%     Y(t0) = [y1(t0), y2(t0), ..., ym(t0)]'.
%
% Inputs:
%     Mn     = matrix M, or function handle for M(t)
%     fcn    = function handle for F(t,Y)
%     Jfcn   = function handle for Jacobian of F, J(t,Y)
%     tvals  = [t0, t1, t2, ..., tN]
%     Y0     = initial value array (column vector of length m)
%     B      = Butcher matrix for IRK coefficients, of the form
%                 B = [c A;
%                      q b;
%                      p b2 ]
%              Here, c is a vector of stage time fractions (s-by-1),
%                    A is a matrix of Butcher coefficients (s-by-s),
%                    q is an integer denoting the method order of accuracy,
%                    b is a vector of solution weights (1-by-s),
%                    p is an integer denoting the embedding order of accuracy,
%                    b2 is a vector of embedding weights (1-by-s),
%              The [p, b2] row is optional.  If that row is not
%              provided the method will default to taking fixed
%              step sizes of size hmin.
%     rtol   = desired relative error of solution  (scalar)
%     atol   = desired absolute error of solution  (vector or scalar)
%     hmin   = minimum internal time step size (hmin <= t(i)-t(i-1), for all i)
%     hmax   = maximum internal time step size (hmax >= hmin)
%     alg    = (optional) solution algorithm:
%                 0 - solve for stages, z_i  [default]
%                 1 - solve for stage RHS, k_i
%
% Outputs: 
%     tvals  = the same as the input array tvals
%     y      = [y(t0), y(t1), y(t2), ..., y(tN)], where each
%               y(t*) is a column vector of length m.
%     nsteps = number of internal time steps taken by method
%     lits   = number of linear solves required by method
%
% Note1: to run in fixed-step mode, call with hmin=hmax as the desired 
% time step size.
%
% Note2: to indicate that the mass matrix _does_not_ depent on time,
% i.e., M ~= M(t), the input "Mn" should be a standard
% double-precision matrix.  To indicate that the mass matrix _does_
% depend on time, the input "Mn" should be a function handle (that
% when called with a scalar-valued "t" returns a standard
% double-precision matrix).
%
% Daniel R. Reynolds
% Department of Mathematics
% Southern Methodist University
% October 2019
% All Rights Reserved


% handle optional 'alg' input
if ~exist('alg','var')
   alg = 0;
end

% set flag based on type of Mn input
MTimeDep = 1;
if (isa(Mn,'double'))
   MTimeDep = 0;
end

% extract Butcher table
[Brows, Bcols] = size(B);
s = Bcols - 1;        % number of stages
c = B(1:s,1);         % stage time fraction array
b = B(s+1,2:s+1)';    % solution weights (convert to column)
A = B(1:s,2:s+1);     % RK coefficients
d = b;                % embedding coefficients (may be overwritten)

% if adaptivity desired, check for embedding coefficients and set the
% order of accuracy accordingly
p = B(Bcols,1);
adaptive = 0;
if (abs(hmax-hmin)/abs(hmax) > sqrt(eps))       % check whether adaptivity is desired
   if (Brows > Bcols)
      if (max(abs(B(Bcols+1,2:Bcols))) > eps)   % check for embedding coeffs
         adaptive = 1;
         p = B(Bcols+1,1);
         d = B(Bcols+1,2:end)';
      end
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
newt_rtol  = rtol/10;      % Newton solver relative tolerance
newt_atol  = atol/10;      % Newton solver absolute tolerance
h_reduce   = 0.1;          % failed step reduction factor 
h_safety   = 0.96;         % adaptivity safety factor
h_growth   = 10;           % adaptivity growth bound
e_bias     = 1.5;          % error bias factor
ONEMSM     = 1-sqrt(eps);  % coefficients to account for
ONEPSM     = 1+sqrt(eps);  %   floating-point roundoff
ERRTOL     = 1.1;          % upper bound on allowed step error
                           %   (in WRMS norm)

% initialize temporary variables
t = tvals(1);
Ynew = Y0;

% create Fdata structure for Newton solver and step solutions
Fdata.MTimeDep = MTimeDep;  % time-depencence of mass matrix
if (MTimeDep)               % mass matrix
   Fdata.M  = Mn(t);
   Fdata.Mn = Mn;
else
   Fdata.M  = Mn;
end
Fdata.f = fcn;              % ODE RHS function handle
Fdata.J = Jfcn;             % ODE RHS Jacobian function handle
Fdata.A = A;                % Butcher table
Fdata.c = c;
Fdata.b = b;
Fdata.d = d;
Fdata.s = s;                % number of stages

% set function names for solve components
if (alg == 1) 
   Init  = @Init_k;   % initializes solution storage
   Guess = @Guess_k;  % initial Newton guess
   Rhs   = @Rhs_k;    % just before equation (45)
   Res   = @Res_k;    % equation (45)
   Jres  = @Jres_k;   % equation (46)
   Sol   = @Sol_k;    % time-evolved solution
   Store = @Store_k;  % stores per-stage results
else
   Init  = @Init_z;   % initializes solution storage
   Guess = @Guess_z;  % initial Newton guess
   Rhs   = @Rhs_z;    % just before equation (41)
   Res   = @Res_z;    % equation (41)
   Jres  = @Jres_z;   % equation (42)
   Sol   = @Sol_z;    % time-evolved solution
   Store = @Store_z;  % stores per-stage results
end

% set functions to compute error weight vector and measure temporal convergence
Ewt = @(Fdata) 1./(rtol*abs(Fdata.yold)+atol);
Rwt = @(Fdata) 1./(rtol*abs(Fdata.M*Fdata.yold)+atol);
WrmsNorm = @(x,w) sqrt(sum((x.*w).^2)/length(x));

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
      if (MTimeDep)      % current mass matrix (M\ used in Init below)
         Fdata.M = Mn(t);
      end

      % set error-weight vectors for this step
      ewt = Ewt(Fdata);
      rwt = Rwt(Fdata);

      % initialize data storage for multiple stages
      [storage,NewtSol] = Init(Y0,Fdata);

      % reset stage failure flag
      st_fail = 0;
      
      % loop over stages
      for stage = 1:s
         
         % Update Fdata structure for current stage
         Fdata.tcur = t + h*c(stage);      % 'time' for current stage
         Fdata.stage = stage;              % current stage index
         Fdata.rhs = Rhs(storage, Fdata);  % 'RHS' of known data
         if (MTimeDep)                     % current mass matrix
            Fdata.M = Mn(Fdata.tcur);
         end
         
         % set nonlinear solver tolerances based on 'alg' type
         if (alg == 1) 
            n_rtol = newt_rtol / h;
            n_atol = newt_atol / h;
         else
            n_rtol = newt_rtol;
            n_atol = newt_atol;
         end         
         
         % set Newton initial guess, call Newton solver, and
         % increment linear solver statistics
         NewtGuess = Guess(NewtSol, Fdata, storage);
         [NewtSol,lin,ierr] = newton(Res, Jres, NewtGuess, Fdata, rwt, newt_maxit);
         lits = lits + lin;
         
         % if Newton method failed, set relevant flags/statistics
         % and break out of stage loop
         if (ierr ~= 0) 
            st_fail = 1;
            c_fails = c_fails + 1;
            break;
         end
         
         % store stage solution
         storage = Store(NewtSol, Fdata, storage);
         
      end
      
      % increment number of internal time steps taken
      nsteps = nsteps + 1;
      
      % compute new solution (and embedding if available)
      [Ynew,Y2] = Sol(storage,Fdata);

      % if stages succeeded and time step adaptivity enabled, check step accuracy
      if ((st_fail == 0) && adaptive)

         % estimate error in current step
	 err_step = e_bias * max(WrmsNorm(Ynew - Y2, ewt), eps);
         
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
         
         % for adaptive methods, use error estimate to adapt the time step
         if (adaptive)

            h_old = h;
            if (err_step == 0.0)     % no error, set max possible
               h = tvals(end)-t;
            else                     % set next h (I-controller)
               h = h_safety * h_old * err_step^(-1.0/p);
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
            error('Cannot achieve desired accuracy.\n  Consider reducing hmin or increasing rtol.\n');
         end

         % otherwise, reset guess, reduce time step, retry solve
         Ynew = Y0;
         h = h * h_reduce;
         
      end  % end logic tests for step success/failure
      
   end  % end while loop attempting to solve steps to next output time

   % store updated solution in output array
   Y(:,tstep) = Ynew;
   
end  % time step loop

% end solve_DIRK_mass function
end



%======= Auxiliary routines when solving for DIRK _stages_ =======%



function [z,NewtGuess] = Init_z(yold, Fdata)
% usage: [z,NewtGuess] = Init_z(yold, Fdata)
%
% Sets aside storage for reusable data in following routines, 
% and initializes first Newton solver guess   
z = zeros(length(yold),Fdata.s);
NewtGuess = yold;
end


function [Y] = Guess_z(Ynew, Fdata, z)
% usage: [Y] = Guess_z(Ynew, Fdata, z)
%
% Sets the initial guess for the Newton iteration stage solve,
% where 'Ynew' was the most-recently-computed solution
Y = Ynew;
end


function [z] = Store_z(Ynew, Fdata, z)
% usage: [z] = Store_z(Ynew, Fdata, z)
%
% Packs reusable data following a stage solve
z(:,Fdata.stage) = Ynew;
end


function [r] = Rhs_z(z, Fdata)
% usage: [r] = Rhs_z(z, Fdata)
%
% Inputs:
%    z     = stage solutions [z_1, ..., z_{stage-1}]
%    Fdata = structure containing extra problem information
%
% Outputs: 
%    r     = rhs vector containing all 'known' information for
%            implicit stage solve
   
% form of r depends on whether M depends on time 
if (Fdata.MTimeDep)
   %    Mi*zi = Mi*y_n + h*Mi*sum_{j=1}^i (a(i,j)*Mj\fj)
   % <=>
   %    Mi*zi - h*(a(i,i)*fi) = Mi*y_n + h*Mi*sum_{j=1}^{i-1} (a(i,j)*Mj\fj)
   % =>
   %    r = Mi*(y_n + h*sum_{j=1}^{i-1} (a(i,j)*Mj\fj))
   r = Fdata.yold;
   for j = 1:Fdata.stage-1
      t = Fdata.t+Fdata.h*Fdata.c(j);
      r = r + Fdata.h*Fdata.A(Fdata.stage,j)*(Fdata.Mn(t) \ Fdata.f(t, z(:,j)));
   end
   r = Fdata.M * r;
else
   %    M*zi = M*y_n + h*sum_{j=1}^i (a(i,j)*fj)
   % <=>
   %    M*zi - h*(a(i,i)*fi) = M*y_n + h*sum_{j=1}^{i-1} (a(i,j)*fj)
   % =>
   %    r = M*y_n + h*sum_{j=1}^{i-1} (a(i,j)*fj)
   r = Fdata.M*Fdata.yold;
   for j = 1:Fdata.stage-1
      r = r + Fdata.h*Fdata.A(Fdata.stage,j)*Fdata.f(Fdata.t+Fdata.h*Fdata.c(j), z(:,j));
   end
end
end


function F = Res_z(z, Fdata)
% usage: F = Res_z(z, Fdata)
%
% Inputs:  z = current guess for stage solution
%          Fdata = structure containing extra information for evaluating F.
% Outputs: F = residual at current guess
%
% This function computes the (non)linear residuals for an intermediate
% stage solution, through calling the user-supplied (in Fdata) ODE
% right-hand side function.
   
F = Fdata.M*z - Fdata.rhs - Fdata.h*Fdata.A(Fdata.stage,Fdata.stage)*Fdata.f(Fdata.tcur, z);

end


function Amat = Jres_z(z, Fdata)
% usage: Amat = Jres_z(z, Fdata)
%
% Inputs:  z = current guess for stage solution
%          Fdata = structure containing extra information for evaluating F.
% Outputs: Amat = Jacobian at current guess
%
% This function computes the Jacobian of each intermediate stage residual
% for a multi-stage DIRK method, through calling the user-supplied (in
% Fdata) ODE Jacobian function. 

Amat = Fdata.M - Fdata.h*Fdata.A(Fdata.stage,Fdata.stage)*Fdata.J(Fdata.tcur, z);

end


function [y,y2] = Sol_z(z, Fdata)
% usage: [y,y2] = Sol_z(z, Fdata)
%
% Inputs:
%    z     = stage solutions [z1, ..., zs]
%    Fdata = structure containing extra problem information
%
% Outputs: 
%    y     = step solution built from the z values
%    y2    = embedded solution (if embedding included in Butcher 
%               table; otherwise the same as y)

% call RHS at each stored stage
f = zeros(size(z,1),Fdata.s);
if (~Fdata.MTimeDep)
   M = Fdata.M;
end
for is=1:Fdata.s
   t = Fdata.t + Fdata.h*Fdata.c(is);
   if (Fdata.MTimeDep)
      M = Fdata.Mn(t);
   end
   f(:,is) = M \ Fdata.f(t, z(:,is));
end

% form the solutions
%    ynew = yold + (h*sum(b(j)*Mj\f(j)))
y  = Fdata.yold + Fdata.h*f*Fdata.b;
y2 = Fdata.yold + Fdata.h*f*Fdata.d;
end



%======= Auxiliary routines when solving for DIRK _RHS_ =======%



function [Kt,NewtGuess] = Init_k(yold, Fdata)
% usage: [Kt,NewtGuess] = Init_k(yold, Fdata)
%
% Sets aside storage for reusable data in following routines, 
% and initializes first Newton solver guess   
Kt = zeros(length(yold),Fdata.s);   % (sol vec) x (stages)
NewtGuess = Fdata.M \ Fdata.f(Fdata.t, yold);
end


function [kt] = Guess_k(ktsol, Fdata, Kt)
% usage: [kt] = Guess_k(ktsol, Fdata, Kt)
%
% Sets the initial guess for the Newton iteration stage solve,
% where 'ksol' was the most-recently-computed solution
kt = ktsol;
end


function [Kt] = Store_k(ktsol, Fdata, Kt)
% usage: [Kt] = Store_k(ktsol, Fdata, Kt)
%
% Packs reusable data following a stage solve
Kt(:,Fdata.stage) = ktsol;
end


function [r] = Rhs_k(Kt, Fdata)
% usage: [r] = Rhs_k(Kt, Fdata)
%
% Inputs:
%    Kt    = stage rhs [M\f(z_1), ..., M\f(z_{stage-1})]
%    Fdata = structure containing extra problem information
%
% Outputs: 
%    r     = rhs vector containing all 'known' information for
%            implicit stage solve
%
%          = y_n + h*sum_{j=1}^{i-1} (a(i,j)*M\fj)
   
% construct rhs
r = Fdata.yold;
for j = 1:Fdata.stage-1
   r = r + Fdata.h*Fdata.A(Fdata.stage,j)*Kt(:,j);
end
end


function F = Res_k(kt, Fdata)
% usage: F = Res_k(kt, Fdata)
%
% Inputs:  kt = current guess for stage rhs
%          Fdata = structure containing extra information for evaluating F.
% Outputs: F = residual at current guess
%
% This function computes the (non)linear residuals for an intermediate
% stage solution, through calling the user-supplied (in Fdata) ODE
% right-hand side function.
   
F = Fdata.M*kt - Fdata.f(Fdata.tcur, Fdata.rhs + Fdata.h*Fdata.A(Fdata.stage,Fdata.stage)*kt);
end


function Amat = Jres_k(kt, Fdata)
% usage: Amat = Jres_k(kt, Fdata)
%
% Inputs:  k = current guess for stage rhs
%          Fdata = structure containing extra information for evaluating F.
% Outputs: Amat = Jacobian at current guess
%
% This function computes the Jacobian of each intermediate stage residual
% for a multi-stage DIRK method, through calling the user-supplied (in
% Fdata) ODE Jacobian function. 

Aii = Fdata.A(Fdata.stage,Fdata.stage);
Amat = Fdata.M - Fdata.h*Aii * Fdata.J(Fdata.tcur, Fdata.rhs + Fdata.h*Aii*kt);
end


function [y,y2] = Sol_k(Kt, Fdata)
% usage: [y,y2] = Sol_k(Kt, Fdata)
%
% Inputs:
%    Kt    = stage RHS vectors [M\f(z1), ..., M\f(zs)]
%    Fdata = structure containing extra problem information
%
% Outputs: 
%    y     = step solution
%    y2    = embedded solution (if embedding included in Butcher 
%               table; otherwise the same as y)

% have RHS at each stored stage, so just piece together
%    ynew = yold + h*sum(b(j)*(Mj\f(j)))
y  = Fdata.yold + Fdata.h*Kt*Fdata.b;
y2 = Fdata.yold + Fdata.h*Kt*Fdata.d;
end
