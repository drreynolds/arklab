function [tvals,Y,nsteps,lits,cfails,afails,ierr] = solve_ARK_bdry(fe,fi,Ji,bdry,tvals,Y0,Be,Bi,rtol,atol,hmin,hmax)
% usage: [tvals,Y,nsteps,lits,cfails,afails,ierr] = solve_ARK_bdry(fe,fi,Ji,bdry,tvals,Y0,Be,Bi,rtol,atol,hmin,hmax)
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
%     bdry   = function to enforce boundary conditions, u = bdry(t,u)
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
%     cfails = number of nonlinear solver convergence failures
%     afails = number of temporal accuracy error failures
%     ierr   = flag denoting success (0) or failure (1)
%
% Note: to run in fixed-step mode, call with hmin=hmax as the desired 
% time step size.
%
% Daniel R. Reynolds
% Department of Mathematics
% Southern Methodist University
% October 2019
% All Rights Reserved

% check for compatible Be,Bi tables
if (size(Be,2) ~= size(Bi,2))   
   error('solve_ARK_bdry error: Be and Bi must have the same number of stages')
end
s = size(Be,2) - 1;          % number of stages
if (Be(s+1,1) ~= Bi(s+1,1))
   error('solve_ARK_bdry error: Be and Bi must have the same method order')
end
if (size(Be,1) > size(Be,2))
   if (Be(s+1,2) ~= Bi(s+1,2))
      error('solve_ARK_bdry error: Be and Bi must have the same embedding order')
   end
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

% initialize outputs
N = length(tvals);
m = length(Y0);
Y = zeros(m,N);
Y(:,1) = Y0;
ierr = 0;

% initialize diagnostics
cfails = 0;   % total convergence failures
afails = 0;   % total accuracy failures

% set the solver parameters
if (adaptive)
   newt_maxit = 3;         % max number of Newton iterations
else
   newt_maxit = 10;        % max number of Newton iterations
end
newt_tol   = 0.1;          % Newton solver tolerance factor
h_cfail    = 0.25;         % failed newton solve step reduction factor 
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

% set functions to compute error weight vector and measure temporal convergence
Ewt = @(Fdata) 1./(rtol*abs(Fdata.yold)+atol);
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

      % set error-weight vector for this step
      ewt = Ewt(Fdata);

      % initialize data storage for multiple stages
      [storage,Fdata] = Init(Y0,Fdata);
      NewtSol = zeros(size(Y0));   % initialize 'correction' solution

      % reset stage failure flag
      st_fail = 0;
      
      % loop over stages
      for stage = 1:s
         
         % update Fdata, set Newton initial guess, and enforce
         % boundary conditions
         Fdata.tcur = t + h*ci(stage);      % 'time' for current [implicit] stage
         Fdata.stage = stage;               % current stage index
         [NewtGuess,Fdata] = Guess(NewtSol, Fdata, bdry, storage);
         Fdata.rhs = Rhs(storage, Fdata);   % 'RHS' of known data
         
         % call Newton solver and increment linear solver statistics
         [NewtSol,lin,nierr] = newton(Res, Jres, NewtGuess, Fdata, ewt, newt_tol, newt_maxit, 0);
         lits = lits + lin;
         
         % if Newton method failed, set relevant flags/statistics
         % and break out of stage loop
         if (nierr ~= 0) 
            st_fail = 1;
            cfails = cfails + 1;
            break;
         end
         
         % enforce boundary conditions and store stage solution
         storage = Store(NewtSol, Fdata, storage);
         
      end
      
      % increment number of internal time steps taken
      nsteps = nsteps + 1;
      
      % compute new solution (and embedding if available), and
      % enforce boundary conditions
      [Ynew,Y2] = Sol(storage,Fdata,bdry);

      % if a stage solve failed
      if (st_fail == 1)

         % if time step adaptivity enabled
         if (adaptive)

            % if already at minimum step, just return with failure
            if (h <= hmin) 
               ierr = 1;
               fprintf('Stage solve failure at minimum step size (t=%g).\n  Consider reducing hmin.\n',Fdata.tcur);
               return
            end

            % otherwise, reset guess, reduce time step, retry solve
            Ynew = Y0;
            h = h * h_cfail;
            continue;
         
         % if time step adaptivity disabled, just return with failure
         else
            ierr = 1;
            fprintf('Stage solve failure in fixed-step mode (t=%g).\n  Consider reducing h.\n',Fdata.tcur);
            return
         end

      end
      
      % if we made it to this point, then all stage solves succeeded
      
      % if time step adaptivity enabled, check step accuracy
      if (adaptive)

         % estimate error in current step
         err_step = e_bias * max(WrmsNorm(Ynew - Y2, ewt), eps);
         
         % if error too high, flag step as a failure (will be be recomputed)
         if (err_step > ERRTOL*ONEPSM) 
            afails = afails + 1;
            st_fail = 1;
            
            % if already at minimum step, just return with failure
            if (h <= hmin) 
               ierr = 1;
               fprintf('Temporal error failure at minimum step size (t=%g).\n  Consider reducing hmin or increasing rtol.\n',Fdata.tcur);
               return
            end
            
         end
         
      end

      % if step was successful (solves succeeded, and error acceptable)
      if (st_fail == 0) 
         
         % update solution and time for last successful step
         Y0 = Ynew;
         t  = t + h;
         
         % for adaptive methods, use error estimate to adapt the time step
         if (adaptive)

            % compute new step size (I-controller)
            h_old = h;
            h = h_safety * h_old * err_step^(-1.0/p);

            % enforce maximum growth rate on step sizes
            h = min(h_growth*h_old, h);

         % otherwise, just use the fixed minimum input step size
         else
            h = hmin;
         end
         
      % if the error test failed
      else

         % if already at minimum step, just return with failure
         if (h <= hmin) 
            ierr = 1;
            fprintf('Cannot achieve desired accuracy at minimum step size (t=%g).\n  Consider reducing hmin or increasing rtol.\n',Fdata.tcur);
            return
         end

         % otherwise, reset guess, reduce time step, retry solve
         Ynew = Y0;
         h_old = h;
         h = min(h_safety * h_old * err_step^(-1.0/p), h_old*h_reduce);
         
      end  % end logic tests for step success/failure
      
   end  % end while loop attempting to solve steps to next output time

   % store updated solution in output array
   Y(:,tstep) = Ynew;
   
end  % time step loop

% end solve_ARK_bdry function
end



%======= Auxiliary routine =======%



function [Z,Fdata] = Init(yold, Fdata)
% usage: [Z,Fdata] = Init(yold, Fdata)
%
% Sets aside storage for stage solutions
Z = zeros(length(yold),Fdata.s);
Fdata.zpred = yold;
end


function [zcor,Fdata] = Guess(zcor, Fdata, bdry, Z)
% usage: [zcor,Fdata] = Guess(zcor, Fdata, bdry, Z)
%
% Sets the initial guess for the Newton iteration stage solve into
% Fdata, where zcor was the most-recently-computed stage correction
Fdata.zpred = Fdata.zpred + zcor;            % update current stored predictor
Fdata.zpred = bdry(Fdata.tcur,Fdata.zpred);  % enforce boundary conditions
zcor = zeros(size(zcor));                    % reset 'correction' guess
end


function [Z] = Store(zcor, Fdata, Z)
% usage: [Z] = Store(zcor, Fdata, Z)
%
% Packs reusable data following a stage solve
Z(:,Fdata.stage) = Fdata.zpred + zcor;
end


function [r] = Rhs(Z, Fdata)
% usage: [r] = Rhs(Z, Fdata)
%
% Inputs:
%    Z     = stage solutions [z_1, ..., z_{stage-1}]
%    Fdata = structure containing extra problem information
%
% Outputs: 
%    r     = rhs vector containing all 'known' information for
%            implicit stage solve
%
%    zi = y_n + h*sum_{j=1}^i (aI(i,j)*fI_j) + h*sum_{j=1}^{i-1} (aE(i,j)*fE_j)
% <=>
%    zi - h*(a(i,i)*fi) = y_n + h*sum_{j=1}^{i-1} (aI(i,j)*fI_j + aE(i,j)*fE_j)
% =>
%    rhs = y_n + h*sum_{j=1}^{i-1} (aI(i,j)*fI_j + aE(i,j)*fE_j)

% construct rhs
r = Fdata.yold;
for j = 1:Fdata.stage-1
   r = r + Fdata.h*Fdata.Ae(Fdata.stage,j)*Fdata.fe(Fdata.t+Fdata.h*Fdata.ce(j), Z(:,j)) ...
         + Fdata.h*Fdata.Ai(Fdata.stage,j)*Fdata.fi(Fdata.t+Fdata.h*Fdata.ci(j), Z(:,j));
end
end


function F = Res(zcor, Fdata)
% usage: F = Res(zcor, Fdata)
%
% Inputs:  zcor = current guess for stage solution correction
%          Fdata = structure containing extra information for evaluating F.
% Outputs: F = residual at current guess
%
% This function computes the (non)linear residuals for an intermediate
% stage solution, through calling the user-supplied (in Fdata) ODE
% right-hand side function.
   
z = Fdata.zpred + zcor;   
F = z - Fdata.rhs - Fdata.h*Fdata.Ai(Fdata.stage,Fdata.stage)*Fdata.fi(Fdata.tcur, z);
end


function Amat = Jres(zcor, Fdata)
% usage: Amat = Jres(zcor, Fdata)
%
% Inputs:  zcor = current guess for stage solution correction
%          Fdata = structure containing extra information for evaluating F.
% Outputs: Amat = Jacobian at current guess
%
% This function computes the Jacobian of each intermediate stage residual
% for a multi-stage ARK method, through calling the user-supplied (in
% Fdata) ODE Jacobian function. 

z = Fdata.zpred + zcor;   
Amat = eye(length(z)) - Fdata.h*Fdata.Ai(Fdata.stage,Fdata.stage)*Fdata.Ji(Fdata.tcur, z);
end


function [y,y2] = Sol(Z, Fdata, bdry)
% usage: [y,y2] = Sol(Z, Fdata, bdry)
%
% Inputs:
%    Z     = stage solutions [z1, ..., zs]
%    Fdata = structure containing extra problem information
%    bdry  = boundary condition enforcement function
%
% Outputs: 
%    y     = step solution built from the Z values
%    y2    = embedded solution (if embedding included in Butcher 
%               table; otherwise the same as y)

% call RHS at each stored stage
fe = zeros(size(Z,1),Fdata.s);
fi = zeros(size(Z,1),Fdata.s);
for is=1:Fdata.s
   t = Fdata.t + Fdata.h*Fdata.ce(is);
   fe(:,is) = Fdata.fe(t, Z(:,is));

   t = Fdata.t + Fdata.h*Fdata.ci(is);
   fi(:,is) = Fdata.fi(t, Z(:,is));
end

% form the solutions
%    ynew = yold + h*sum(be(j)*fe(j) + bi(j)*fi(j))
y  = Fdata.yold + Fdata.h*fe*Fdata.be + Fdata.h*fi*Fdata.bi;
y2 = Fdata.yold + Fdata.h*fe*Fdata.de + Fdata.h*fi*Fdata.di;

% enforce boundary conditions on both solutions
y  = bdry(Fdata.t+Fdata.h,y);
y2 = bdry(Fdata.t+Fdata.h,y2);

end
