function [tvals,Y,nsteps,lits] = solve_DIRK(fcn,Jfcn,tvals,Y0,B,rtol,atol,hmin,hmax,alg)
% usage: [tvals,Y,nsteps,lits] = solve_DIRK(fcn,Jfcn,tvals,Y0,B,rtol,atol,hmin,hmax,alg)
%
% Adaptive time step diagonally-implicit Runge-Kutta solver for the
% vector-valued ODE problem 
%     y' = F(t,Y), t in tvals, y in R^m,
%     Y(t0) = [y1(t0), y2(t0), ..., ym(t0)]'.
%
% Inputs:
%     fcn    = function handle for F(t,Y)
%     Jfcn   = function handle for Jacobian of F, J(t,Y)%
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
% Note: to run in fixed-step mode, call with hmin=hmax as the desired 
% time step size, and set the tolerances to large positive numbers.
%
% Daniel R. Reynolds
% Department of Mathematics
% Southern Methodist University
% March 2017
% All Rights Reserved


% handle optional 'alg' input
if ~exist('alg','var')
   alg = 0;
end

% extract DIRK method information from B
[Brows, Bcols] = size(B);
s = Bcols - 1;        % number of stages
c = B(1:s,1);         % stage time fraction array
b = (B(s+1,2:s+1))';  % solution weights (convert to column)
A = B(1:s,2:s+1);     % RK coefficients
q = B(s+1,1);         % method order

% initialize as non-embedded, until proven otherwise
embedded = 0;
p = 0;
if (Brows > Bcols)
   if (max(abs(B(s+2,2:s+1))) > eps)
      embedded = 1;
      d = (B(s+2,2:s+1))';
      p = B(s+2,1);
   end
else
   d = b;
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
Fdata.f = fcn;    % ODE RHS function handle
Fdata.J = Jfcn;   % ODE RHS Jacobian function handle
Fdata.A = A;      % Butcher table
Fdata.c = c;
Fdata.b = b;
Fdata.d = d;
Fdata.s = s;      % number of stages

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
      [storage,NewtSol] = Init(Y0,Fdata);

      % reset stage failure flag
      st_fail = 0;
      
      % loop over stages
      for stage = 1:s
         
         % Update Fdata structure for current stage
         Fdata.tcur = t + h*c(stage);      % 'time' for current stage
         Fdata.stage = stage;              % current stage index
         Fdata.rhs = Rhs(storage, Fdata);  % 'RHS' of known data

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
         [NewtSol,lin,ierr] = newton(Res, Jres, NewtGuess, Fdata, n_rtol, n_atol, newt_maxit);
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

% end solve_DIRK function
end



%======= Auxiliary routines when solving for ARK _stages_ =======%



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
%
%    zi = y_n + h*sum_{j=1}^i (a(i,j)*fj)
% <=>
%    zi - h*(a(i,i)*fi) = y_n + h*sum_{j=1}^{i-1} (a(i,j)*fj)
% =>
%    rhs = y_n + h*sum_{j=1}^{i-1} (a(i,j)*fj)
   
% construct rhs
r = Fdata.yold;
for j = 1:Fdata.stage-1
   r = r + Fdata.h*Fdata.A(Fdata.stage,j)*Fdata.f(Fdata.t+Fdata.h*Fdata.c(j), z(:,j));
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
   
F = z - Fdata.rhs - Fdata.h*Fdata.A(Fdata.stage,Fdata.stage)*Fdata.f(Fdata.tcur, z);

end


function Amat = Jres_z(z, Fdata)
% usage: Amat = Jres_z(z, Fdata)
%
% Inputs:  z = current guess for stage solution
%          Fdata = structure containing extra information for evaluating F.
% Outputs: Amat = Jacobian at current guess
%
% This function computes the Jacobian of each intermediate stage residual
% for a multi-stage ARK method, through calling the user-supplied (in
% Fdata) ODE Jacobian function. 

Amat = eye(length(z)) - Fdata.h*Fdata.A(Fdata.stage,Fdata.stage)*Fdata.J(Fdata.tcur, z);

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
for is=1:Fdata.s
   t = Fdata.t + Fdata.h*Fdata.c(is);
   f(:,is) = Fdata.f(t, z(:,is));
end

% form the solutions
%    ynew = yold + h*sum(b(j)*f(j))
y  = Fdata.yold + Fdata.h*f*Fdata.b;
y2 = Fdata.yold + Fdata.h*f*Fdata.d;
end



%======= Auxiliary routines when solving for ARK _RHS_ =======%



function [K,NewtGuess] = Init_k(yold, Fdata)
% usage: [K,NewtGuess] = Init_k(yold, Fdata)
%
% Sets aside storage for reusable data in following routines, 
% and initializes first Newton solver guess   
K = zeros(length(yold),Fdata.s);   % (sol vec) x (stages)
NewtGuess = Fdata.f(Fdata.t, yold);
end


function [k] = Guess_k(ksol, Fdata, K)
% usage: [k] = Guess_k(ksol, Fdata, K)
%
% Sets the initial guess for the Newton iteration stage solve,
% where 'ksol' was the most-recently-computed solution
k = ksol;
end


function [K] = Store_k(ksol, Fdata, K)
% usage: [K] = Store_k(ksol, Fdata, K)
%
% Packs reusable data following a stage solve
K(:,Fdata.stage) = ksol;
end


function [r] = Rhs_k(K, Fdata)
% usage: [r] = Rhs_k(K, Fdata)
%
% Inputs:
%    K     = stage rhs [f(z_1), ..., f(z_{stage-1})]
%    Fdata = structure containing extra problem information
%
% Outputs: 
%    r     = rhs vector containing all 'known' information for
%            implicit stage solve
%
%    rhs = y_n + h*sum_{j=1}^{i-1} (a(i,j)*fj)
   
% construct rhs
r = Fdata.yold;
for j = 1:Fdata.stage-1
   r = r + Fdata.h*Fdata.A(Fdata.stage,j)*K(:,j);
end
end


function F = Res_k(k, Fdata)
% usage: F = Res_k(k, Fdata)
%
% Inputs:  k = current guess for stage rhs
%          Fdata = structure containing extra information for evaluating F.
% Outputs: F = residual at current guess
%
% This function computes the (non)linear residuals for an intermediate
% stage solution, through calling the user-supplied (in Fdata) ODE
% right-hand side function.
   
F = k - Fdata.f(Fdata.tcur, Fdata.rhs + Fdata.h*Fdata.A(Fdata.stage,Fdata.stage)*k);
end


function Amat = Jres_k(k, Fdata)
% usage: Amat = Jres_k(k, Fdata)
%
% Inputs:  k = current guess for stage rhs
%          Fdata = structure containing extra information for evaluating F.
% Outputs: Amat = Jacobian at current guess
%
% This function computes the Jacobian of each intermediate stage residual
% for a multi-stage ARK method, through calling the user-supplied (in
% Fdata) ODE Jacobian function. 

Aii = Fdata.A(Fdata.stage,Fdata.stage);
Amat = eye(length(k)) - Fdata.h*Aii * Fdata.J(Fdata.tcur, Fdata.rhs + Fdata.h*Aii*k);
end


function [y,y2] = Sol_k(K, Fdata)
% usage: [y,y2] = Sol_k(K, Fdata)
%
% Inputs:
%    K     = stage RHS vectors [f(z1), ..., f(zs)]
%    Fdata = structure containing extra problem information
%
% Outputs: 
%    y     = step solution
%    y2    = embedded solution (if embedding included in Butcher 
%               table; otherwise the same as y)

% have RHS at each stored stage, so just piece together
%    ynew = yold + h*sum(b(j)*f(j))
y  = Fdata.yold + Fdata.h*K*Fdata.b;
y2 = Fdata.yold + Fdata.h*K*Fdata.d;
end
