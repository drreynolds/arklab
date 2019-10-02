function [tvals, Y, nsteps, lits] = solve_ARK_mass(Mn, fe, fi, Ji, tvals, Y0, Be, ...
                                                   Bi, rtol, atol, hmin, hmax, alg)
% usage: [tvals, Y, nsteps, lits] = solve_ARK_mass(Mn, fe, fi, Ji, tvals, Y0, Be, ...
%                                                  Bi, rtol, atol, hmin, hmax, alg)
%
% Adaptive time step additive Runge-Kutta solver for the
% vector-valued ODE problem  
%     M(t) * y'(t) = fe(t,Y) + fi(t,Y), t in tvals, y in R^m,
%     Y(t0) = [y1(t0), y2(t0), ..., ym(t0)]'.
%
% Inputs:
%     Mn     = matrix M, or function handle for M(t)
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

% check for compatible Be,Bi tables
if (size(Be,2) ~= size(Bi,2))   
   error('solve_ARK_mass error: Be and Bi must have the same number of stages')
end
s = size(Be,2) - 1;          % number of stages
if (Be(s+1,1) ~= Bi(s+1,1))
   error('solve_ARK_mass error: Be and Bi must have the same method order')
end
if (size(Be,1) > size(Be,2))
   if (Be(s+1,2) ~= Bi(s+1,2))
      error('solve_ARK_mass error: Be and Bi must have the same embedding order')
   end
end

% set flag based on type of Mn input
MTimeDep = 1;
if (isa(Mn,'double'))
   MTimeDep = 0;
end

% extract ERK method information from Be
[Brows, Bcols] = size(Be);
ce = Be(1:s,1);         % stage time fraction array
be = Be(s+1,2:s+1)';    % solution weights (convert to column)
Ae = Be(1:s,2:s+1);     % RK coefficients
de = be;                % embedding coefficients (may be overwritten)

% extract DIRK method information from Bi
[Brows, Bcols] = size(Bi);
ci = Bi(1:s,1);         % stage time fraction array
bi = Bi(s+1,2:s+1)';    % solution weights (convert to column)
Ai = Bi(1:s,2:s+1);     % RK coefficients
di = bi;                % embedding coefficients (may be overwritten)

% if using a time-dependent mass matrix, ensure that ci==ce
if (MTimeDep && (norm(ce-ci,inf)>100*eps))
   error('solve_ARK_mass error: ce and ci must match when using a time-dependent mass matrix')
end

% if adaptivity desired, check for embedding coefficients and set the
% order of accuracy accordingly
p = Be(s+1,1);
adaptive = 0;
if (abs(hmax-hmin)/abs(hmax) > sqrt(eps))       % check whether adaptivity is desired
   if (Brows > Bcols)
      if ( (max(abs(Be(s+2,2:s+1))) > eps) && ...
           (max(abs(Bi(s+2,2:s+1))) > eps) )    % check embedding coeffs
         adaptive = 1;
         p = Be(s+2,1);
         de = Be(s+2,2:s+1)';
         di = Bi(s+2,2:s+1)';
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
Fdata.fe = fe;    % ODE RHS function names
Fdata.fi = fi;
Fdata.Ji = Ji;    % ODE RHS Jacobian function name
Fdata.Ae = Ae;    % Butcher tables
Fdata.ce = ce;
Fdata.be = be;
Fdata.de = de;
Fdata.Ai = Ai;
Fdata.ci = ci;
Fdata.bi = bi;
Fdata.di = di;
Fdata.s  = s;     % number of stages

% set function names for solve components, depending on the choice of 'alg'
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
   while ((t-tvals(tstep))*h < 0)
      
      % bound internal time step 
      h = max([h, hmin]);            % enforce minimum time step size
      h = min([h, hmax]);            % maximum time step size
      h = min([h, tvals(tstep)-t]);  % stop at output time

      % set Fdata values for this step
      Fdata.h    = h;    % current step size
      Fdata.yold = Y0;   % solution from previous step
      Fdata.t    = t;    % time of last successful step
      if (MTimeDep)      % current mass matrix (M^{-1} used in Init below)
         Fdata.M = Mn(t);
      end

      % set error-weight vector for this step
      ewt = Ewt(Fdata);
      rwt = Rwt(Fdata);

      % initialize data storage for multiple stages
      [storage,NewtSol] = Init(Y0,Fdata);

      % reset stage failure flag
      st_fail = 0;
      
      % loop over stages
      for stage = 1:s

         % Update Fdata structure for current stage
         Fdata.tcur = t + h*ci(stage);     % 'time' for current [implicit] stage
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
      if ((st_fail == 0) & adaptive)

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
         
         % for embedded methods, use error estimate to adapt the time step
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

% end solve_ARK_mass function
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

% form of r depends on whether M depends on time 
if (Fdata.MTimeDep)
   %    Mi*zi = Mi*y_n + h*Mi*sum_{j=1}^i Mj^{-1}*(aI(i,j)*fI(j) + aE(i,j)*fE(j))
   % <=>
   %    Mi*zi - h*(aI(i,i)*fI(i)) = Mi*y_n + h*Mi*sum_{j=1}^{i-1} Mj^{-1}*(aI(i,j)*fI(j) + aE(i,j)*fE(j))
   % =>
   %    r = Mi*(y_n + h*( sum_{j=1}^{i-1} Mj^{-1}*(aI(i,j)*fI(j)) + aE(i,j)*fE(j)))
   r = Fdata.yold;
   for j = 1:Fdata.stage-1
      t = Fdata.t+Fdata.h*Fdata.ce(j);  % note that ci=ce is required
      Mj = Fdata.Mn(t);
      r = r + Fdata.h * (Mj \ ( Fdata.Ae(Fdata.stage,j)*Fdata.fe(t, z(:,j)) ...
                              + Fdata.Ai(Fdata.stage,j)*Fdata.fi(t, z(:,j)) ));
   end
   r = Fdata.M * r;
else
   %    M*zi = M*y_n + h*sum_{j=1}^i (aI(i,j)*fI_j) + h*sum_{j=1}^{i-1} (aE(i,j)*fE_j)
   % <=>
   %    M*zi - h*(a(i,i)*fi) = M*y_n + h*sum_{j=1}^{i-1} (aI(i,j)*fI_j + aE(i,j)*fE_j)
   % =>
   %    r = M*y_n + h*sum_{j=1}^{i-1} (aI(i,j)*fI_j + aE(i,j)*fE_j)
   r = Fdata.M*Fdata.yold;
   for j = 1:Fdata.stage-1
      r = r + Fdata.h * ( Fdata.Ae(Fdata.stage,j)*Fdata.fe(Fdata.t+Fdata.h*Fdata.ce(j), z(:,j)) ...
                        + Fdata.Ai(Fdata.stage,j)*Fdata.fi(Fdata.t+Fdata.h*Fdata.ci(j), z(:,j)) );
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
   
F = Fdata.M*z - Fdata.rhs - Fdata.h*Fdata.Ai(Fdata.stage,Fdata.stage)*Fdata.fi(Fdata.tcur, z);

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

Amat = Fdata.M - Fdata.h*Fdata.Ai(Fdata.stage,Fdata.stage)*Fdata.Ji(Fdata.tcur, z);

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
fe = zeros(size(z,1),Fdata.s);
fi = zeros(size(z,1),Fdata.s);
if (~Fdata.MTimeDep)
   M = Fdata.M;
end
for is=1:Fdata.s
   if (Fdata.MTimeDep)
      t = Fdata.t + Fdata.h*Fdata.ce(is);  % note: ce=ci is required
      Mj = Fdata.Mn(t);
      fe(:,is) = Mj \ Fdata.fe(t, z(:,is));
      fi(:,is) = Mj \ Fdata.fi(t, z(:,is));
   else
      t = Fdata.t + Fdata.h*Fdata.ce(is);
      fe(:,is) = Fdata.M \ Fdata.fe(t, z(:,is));
      
      t = Fdata.t + Fdata.h*Fdata.ci(is);
      fi(:,is) = Fdata.M \ Fdata.fi(t, z(:,is));
   end
end

% form the solutions
%    ynew = yold + h*sum(be(j)*Mj\fe(j) + bi(j)*Mj\fi(j))
y  = Fdata.yold + Fdata.h*fe*Fdata.be + Fdata.h*fi*Fdata.bi;
y2 = Fdata.yold + Fdata.h*fe*Fdata.de + Fdata.h*fi*Fdata.di;
end



%======= Auxiliary routines when solving for ARK _RHS_ =======%



function [Kt,NewtGuess] = Init_k(yold, Fdata)
% usage: [Kt,NewtGuess] = Init_k(yold, Fdata)
%
% Sets aside storage for reusable data in following routines, 
% and initializes first Newton solver guess   
Kt = zeros(length(yold),Fdata.s,2);   % (sol vec) x (stages) x (implicit,explicit)
NewtGuess = Fdata.M \ Fdata.fi(Fdata.t, yold);
end


function [kt] = Guess_k(ktsol, Fdata, K)
% usage: [kt] = Guess_k(ktsol, Fdata, K)
%
% Sets the initial guess for the Newton iteration stage solve,
% where 'ksol' was the most-recently-computed solution
kt = ktsol;
end


function [Kt] = Store_k(ktsol, Fdata, Kt)
% usage: [Kt] = Store_k(ktsol, Fdata, Kt)
%
% Packs reusable data following a stage solve

% just store implicit RHS
Kt(:,Fdata.stage,1) = ktsol;

% compute/store explicit RHS at new stage
znew = Fdata.rhs + Fdata.h*Fdata.Ai(Fdata.stage,Fdata.stage)*ktsol;
Kt(:,Fdata.stage,2) = Fdata.M \ Fdata.fe(Fdata.t+Fdata.h*Fdata.ce(Fdata.stage), znew);
end


function [r] = Rhs_k(Kt, Fdata)
% usage: [r] = Rhs_k(Kt, Fdata)
%
% Inputs:
%    Kt    = stage rhs [M1^{-1}*fi(z_1), ..., M_{st-1}^{-1}*fi(z_{st-1}), M1^{-1}*fe(z_1), ..., M_{st-1}^{-1}*fe(z_{st-1})]
%    Fdata = structure containing extra problem information
%
% Outputs: 
%    r     = rhs vector containing all 'known' information for
%            implicit stage solve
%
%    rhs = y_n + h*sum_{j=1}^{i-1} (aI(i,j)*Mj^{-1}*fI_j + aE(i,j)*Mj^{-1}*fE_j)

% construct rhs
r = Fdata.yold;
for j = 1:Fdata.stage-1
   r = r + Fdata.h*(Fdata.Ai(Fdata.stage,j)*Kt(:,j,1) + Fdata.Ae(Fdata.stage,j)*Kt(:,j,2));
end
end


function F = Res_k(kt, Fdata)
% usage: F = Res_k(kt, Fdata)
%
% Inputs:  kt = current guess for implicit stage rhs
%          Fdata = structure containing extra information for evaluating F.
% Outputs: F = residual at current guess
%
% This function computes the (non)linear residuals for an intermediate
% stage solution, through calling the user-supplied (in Fdata) ODE
% right-hand side function.
   
F = Fdata.M*kt - Fdata.fi(Fdata.tcur, Fdata.rhs + Fdata.h*Fdata.Ai(Fdata.stage,Fdata.stage)*kt);
end


function Amat = Jres_k(kt, Fdata)
% usage: Amat = Jres_k(kt, Fdata)
%
% Inputs:  kt = current guess for implicit stage rhs
%          Fdata = structure containing extra information for evaluating F.
% Outputs: Amat = Jacobian at current guess
%
% This function computes the Jacobian of each intermediate stage residual
% for a multi-stage ARK method, through calling the user-supplied (in
% Fdata) ODE Jacobian function. 

Aii = Fdata.Ai(Fdata.stage,Fdata.stage);
Amat = Fdata.M - Fdata.h*Aii*Fdata.Ji(Fdata.tcur, Fdata.rhs + Fdata.h*Aii*kt);
end


function [y,y2] = Sol_k(Kt, Fdata)
% usage: [y,y2] = Sol_k(Kt, Fdata)
%
% Inputs:
%    Kt    = stage rhs [M\fi(z_1), ..., M\fi(z_{stage-1}), M\fe(z_1), ..., M\fe(z_{stage-1})]
%    Fdata = structure containing extra problem information
%
% Outputs: 
%    y     = step solution
%    y2    = embedded solution (if embedding included in Butcher 
%               table; otherwise the same as y)

% have RHS at each stored stage, so just piece together
%    ynew = yold + h*sum(bi(j)*Mj^{-1}*fi(j) + be(j)*Mj^{-1}*fe(j))
y  = Fdata.yold + Fdata.h*Kt(:,:,1)*Fdata.bi + Fdata.h*Kt(:,:,2)*Fdata.be;
y2 = Fdata.yold + Fdata.h*Kt(:,:,1)*Fdata.di + Fdata.h*Kt(:,:,2)*Fdata.de;
end
