% Driver for heat equation with time-dependent Dirichlet boundary conditions:
%      u_t = lambda*u_xx + f,  0<x<pi,  0<t<4/lambda
%      u(0,x) = 0,             0<x<pi,
%      u(t,0) = b1(t),                  0<t<4/lambda
%      u(t,pi) = b2(t),                 0<t<4/lambda
% The forcing and boundary terms have the form
%    f(t,x) = \sum_{k=1}^M c_k \cos(kx)
%    b1(t)  = \sum_{k=1}^M c_k/(lambda*k^2) (1 - e^{-t*lambda*k^2})
%    b2(t)  = \sum_{k=1}^M c_k/(lambda*k^2) (1 - e^{-t*lambda*k^2}) \cos(k\pi)
% and thus u has analytical solution
%    u(t,x) = \sum_{k=1}^M c_k/(lambda*k^2) (1 - e^{-t*lambda*k^2}) \cos(kx).
%
% Since the highest modes equilibrate the fastest, we use
% real-valued coefficients c_k such that
%
%    |c_1| < |c_2| < ... < |c_M|
%
% Daniel R. Reynolds
% Department of Mathematics
% Southern Methodist University
% November 2019
% All Rights Reserved
clear

% general driver flags
create_diary = true;
store_output = true;
run_natural = true;
run_forced = true;
run_approx = true;

% start output diary
if (create_diary)
  !\rm timedep_bdry_results.txt
  diary timedep_bdry_results.txt
end

% set problem parameters
lambdas = [1, 10, 100, 1000, 1e4, 1e5];
xl = 0;
xr = pi;
m = 50;
T0 = 0;
global Pdata;
Pdata.c = [1, 2, 3, 4, 5];
Pdata.m = m;
Pdata.xspan = linspace(xl,xr,m)';
Pdata.dx = (Pdata.xspan(end)-Pdata.xspan(1))/(m-1);

% set testing parameters
hvals = [0.5, 0.25, 0.1, 0.05, 0.025, 0.01, 0.005];
rtol = 1e-7;
atol = 1e-13*ones(m,1);
Y0 = zeros(size(Pdata.xspan));
DIRKmethods = {'SDIRK-2-2','EDIRK-3-3','Kvaerno(5,3,4)-ESDIRK','Kvaerno(7,4,5)-ESDIRK'};
ARKEmethods = {'ARK(2,3,2)-ERK','Ascher(2,3,3)-ERK','ARK4(3)7L[2]SA-ERK','ARK5(4)8L[2]SA-ERK'};
ARKImethods = {'ARK(2,3,2)-SDIRK','Ascher(2,3,3)-SDIRK','ARK4(3)7L[2]SA-ESDIRK','ARK5(4)8L[2]SA-ESDIRK'};
ERKmethods  = {'ERK-1-1','ERK-2-2','Ascher(2,3,3)-ERK','ERK-4-4','Dormand-Prince-ERK'};


% general output header information
fprintf('\nTime-dependent boundary test problem 2:\n')
fprintf('    spatial domain: [%g, %g]\n',Pdata.xspan(1),Pdata.xspan(end))
fprintf('    spatial mesh nodes: %i\n',m)
fprintf('    tolerances:  rtol = %g,  atol = %g\n', rtol, atol(1));


% create 'statistics' storage
DIRK_natural_rmserrs = zeros(length(lambdas),length(DIRKmethods),length(hvals));
DIRK_natural_ord = zeros(length(lambdas),length(DIRKmethods));
DIRK_natural_ordred = zeros(length(lambdas),length(DIRKmethods));
DIRK_forced_rmserrs = zeros(length(lambdas),length(DIRKmethods),length(hvals));
DIRK_forced_ord = zeros(length(lambdas),length(DIRKmethods));
DIRK_forced_ordred = zeros(length(lambdas),length(DIRKmethods));
DIRK_approx_rmserrs = zeros(length(lambdas),length(DIRKmethods),length(hvals));
DIRK_approx_ord = zeros(length(lambdas),length(DIRKmethods));
DIRK_approx_ordred = zeros(length(lambdas),length(DIRKmethods));

ARK_natural_rmserrs = zeros(length(lambdas),length(ARKEmethods),length(hvals));
ARK_natural_ord = zeros(length(lambdas),length(ARKEmethods));
ARK_natural_ordred = zeros(length(lambdas),length(ARKEmethods));
ARK_forced_rmserrs = zeros(length(lambdas),length(ARKEmethods),length(hvals));
ARK_forced_ord = zeros(length(lambdas),length(ARKEmethods));
ARK_forced_ordred = zeros(length(lambdas),length(ARKEmethods));
ARK_approx_rmserrs = zeros(length(lambdas),length(ARKEmethods),length(hvals));
ARK_approx_ord = zeros(length(lambdas),length(ARKEmethods));
ARK_approx_ordred = zeros(length(lambdas),length(ARKEmethods));

ERK_natural_rmserrs = zeros(length(lambdas),length(ERKmethods),length(hvals));
ERK_natural_ord = zeros(length(lambdas),length(ERKmethods));
ERK_natural_ordred = zeros(length(lambdas),length(ERKmethods));
ERK_forced_rmserrs = zeros(length(lambdas),length(ERKmethods),length(hvals));
ERK_forced_ord = zeros(length(lambdas),length(ERKmethods));
ERK_forced_ordred = zeros(length(lambdas),length(ERKmethods));
ERK_approx_rmserrs = zeros(length(lambdas),length(ERKmethods),length(hvals));
ERK_approx_ord = zeros(length(lambdas),length(ERKmethods));
ERK_approx_ordred = zeros(length(lambdas),length(ERKmethods));

% loop over stiffness parameters
for il = 1:length(lambdas)

   % set lambda-specific testing values and output general testing information
   lambda = lambdas(il);
   Tf = 5/lambda;
   tout = linspace(T0,Tf,11);
   Pdata.lambda = lambda;
   fprintf('\nStiffness factor lambda = %g\n', lambda);

   % set problem-defining functions  (t must be scalar in each)
   fn = @f_timedep_bdry;
   fe = @fe_timedep_bdry;
   fi = @fi_timedep_bdry;
   Jn = @J_timedep_bdry;
   Ji = @J_timedep_bdry;
   Es = @(t,y) Tf;
   k = 1:length(Pdata.c);
   Pdata.b1 = @(t) ( cos(xl*1)*(Pdata.c(1))/(Pdata.lambda*1)*(1-exp(-Pdata.lambda*t)) ...
                   + cos(xl*2)*(Pdata.c(2))/(Pdata.lambda*4)*(1-exp(-Pdata.lambda*t*4)) ...
                   + cos(xl*3)*(Pdata.c(3))/(Pdata.lambda*9)*(1-exp(-Pdata.lambda*t*9)) ...
                   + cos(xl*4)*(Pdata.c(4))/(Pdata.lambda*16)*(1-exp(-Pdata.lambda*t*16)) ...
                   + cos(xl*5)*(Pdata.c(5))/(Pdata.lambda*25)*(1-exp(-Pdata.lambda*t*25)) );
   Pdata.b2 = @(t) ( cos(xr*1)*(Pdata.c(1))/(Pdata.lambda*1)*(1-exp(-Pdata.lambda*t)) ...
                   + cos(xr*2)*(Pdata.c(2))/(Pdata.lambda*4)*(1-exp(-Pdata.lambda*t*4)) ...
                   + cos(xr*3)*(Pdata.c(3))/(Pdata.lambda*9)*(1-exp(-Pdata.lambda*t*9)) ...
                   + cos(xr*4)*(Pdata.c(4))/(Pdata.lambda*16)*(1-exp(-Pdata.lambda*t*16)) ...
                   + cos(xr*5)*(Pdata.c(5))/(Pdata.lambda*25)*(1-exp(-Pdata.lambda*t*25)) );
   Pdata.b1dot = @(t) ( cos(xl*1)*(Pdata.c(1))*(exp(-Pdata.lambda*t)) ...
                      + cos(xl*2)*(Pdata.c(2))*(exp(-Pdata.lambda*t*4)) ...
                      + cos(xl*3)*(Pdata.c(3))*(exp(-Pdata.lambda*t*9)) ...
                      + cos(xl*4)*(Pdata.c(4))*(exp(-Pdata.lambda*t*16)) ...
                      + cos(xl*5)*(Pdata.c(5))*(exp(-Pdata.lambda*t*25)) );
   Pdata.b2dot = @(t) ( cos(xr*1)*(Pdata.c(1))*(exp(-Pdata.lambda*t)) ...
                      + cos(xr*2)*(Pdata.c(2))*(exp(-Pdata.lambda*t*4)) ...
                      + cos(xr*3)*(Pdata.c(3))*(exp(-Pdata.lambda*t*9)) ...
                      + cos(xr*4)*(Pdata.c(4))*(exp(-Pdata.lambda*t*16)) ...
                      + cos(xr*5)*(Pdata.c(5))*(exp(-Pdata.lambda*t*25)) );

   % generate reference solution
   fprintf('\nGenerating reference solution with ode15s:\n')
   opts = odeset('RelTol',1e-13, 'AbsTol',atol,'Jacobian',Jn);
   [t,Yref] = ode15s(fn, tout, Y0, opts);
   Yref = Yref';

   % solve with 'natural' approach
   if (run_natural)
     fprintf('\n  Results with natural approach:\n')

     % run DIRK tests
     if (length(DIRKmethods) > 0)
       for ib = 1:length(DIRKmethods)

         mname = DIRKmethods{ib};
         B = butcher(mname);  s = numel(B(1,:))-1;  q = B(s+1,1);
         fprintf('\n    DIRK integrator: %s (order = %i)\n',mname,q)
         errs = zeros(size(hvals));
         order = zeros(length(hvals)-1,1);
         for ih = 1:length(hvals)
           hval = hvals(ih)/lambda;
           fprintf('      h = %.5e,',hval);
           [t,Y,ns,nl,cf,af] = solve_DIRK(fn, Jn, tout, Y0, B, rtol, atol, hval, hval, 0);
           errs(ih) = sqrt(sum(sum((Y-Yref).^2))/numel(Y));
           if (isnan(errs(ih)) || isinf(errs(ih)) || errs(ih) > 10)
             errs(ih) = 1;
           end
           DIRK_natural_rmserrs(il,ib,ih) = errs(ih);
           if (ih>1)
             order(ih-1) = log( errs(ih)/errs(ih-1) ) / log( hvals(ih)/hvals(ih-1) );
             fprintf('  rmserr = %.2e,  lsolves = %i,  order = %.2e\n',errs(ih),nl,order(ih-1));
           else
             fprintf('  rmserr = %.2e,  lsolves = %i\n',errs(ih),nl);
           end
           fprintf('      newton conv fails = %i, temporal error fails = %i\n',cf,af);
         end
         s = sort(order);  ord = s(end-1);
         fprintf('    estimated order = %g,  reduction = %g\n',ord,max(0,q-ord))
         DIRK_natural_ord(il,ib) = ord;
         DIRK_natural_ordred(il,ib) = max(0,q-ord);

       end
     end

     % run ARK tests
     if ((length(ARKImethods) == length(ARKEmethods)) && (length(ARKImethods)>0))
       for ib = 1:length(ARKImethods)

         mname1 = ARKEmethods{ib};
         Be = butcher(mname1);  s = numel(Be(1,:))-1;  q = Be(s+1,1);
         mname2 = ARKImethods{ib};
         Bi = butcher(mname2);
         fprintf('\n    ARK integrator: %s/%s (order = %i)\n',mname1,mname2,q)
         errs = zeros(size(hvals));
         order = zeros(length(hvals)-1,1);
         for ih = 1:length(hvals)
           hval = hvals(ih)/lambda;
           fprintf('      h = %.5e,',hval);
           [t,Y,ns,nl,cf,af] = solve_ARK(fe, fi, Ji, tout, Y0, Be, Bi, rtol, atol, hval, hval, 0);
           errs(ih) = sqrt(sum(sum((Y-Yref).^2))/numel(Y));
           if (isnan(errs(ih)) || isinf(errs(ih)) || errs(ih) > 10)
             errs(ih) = 1;
           end
           ARK_natural_rmserrs(il,ib,ih) = errs(ih);
           if (ih>1)
             order(ih-1) = log( errs(ih)/errs(ih-1) ) / log( hvals(ih)/hvals(ih-1) );
             fprintf('  rmserr = %.2e,  lsolves = %i,  order = %.2e\n',errs(ih),nl,order(ih-1));
           else
             fprintf('  rmserr = %.2e,  lsolves = %i\n',errs(ih),nl);
           end
           fprintf('      newton conv fails = %i, temporal error fails = %i\n',cf,af);
         end
         s = sort(order);  ord = s(end-1);
         fprintf('    estimated order = %g,  reduction = %g\n',ord,max(0,q-ord))
         ARK_natural_ord(il,ib) = ord;
         ARK_natural_ordred(il,ib) = max(0,q-ord);

       end
     end

     % run ERK tests
     if ((length(ERKmethods) > 0) && (abs(lambda)<=100))
       for ib = 1:length(ERKmethods)
         mname = ERKmethods{ib};
         B = butcher(mname);  s = numel(B(1,:))-1;  q = B(s+1,1);
         fprintf('\n   ERK integrator: %s (order = %i)\n',mname,q)
         errs_rms = zeros(size(hvals));
         errs_max = zeros(size(hvals));
         order = zeros(length(hvals)-1,1);
         for ih = 1:length(hvals)
           hval = hvals(ih)/lambda/20;
           fprintf('      h = %.5e,',hval);
           [t,Y,ns,af] = solve_ERK(fn, Es, tout, Y0, B, rtol, atol, hval, hval);
           errs(ih) = sqrt(sum(sum((Y-Yref).^2))/numel(Y));
           if (isnan(errs(ih)) || isinf(errs(ih)) || errs(ih) > 10)
             errs(ih) = 1;
           end
           ARK_natural_rmserrs(il,ib,ih) = errs(ih);
           if (ih>1)
             order(ih-1) = log( errs(ih)/errs(ih-1) ) / log( hvals(ih)/hvals(ih-1) );
             fprintf('  rmserr = %.2e,  lsolves = %i,  order = %.2e\n',errs(ih),nl,order(ih-1));
           else
             fprintf('  rmserr = %.2e,  lsolves = %i\n',errs(ih),nl);
           end
           fprintf('      temporal error fails = %i\n',af);
         end
         s = sort(order);  ord = s(end-1);
         fprintf('    estimated order = %g,  reduction = %g\n',ord,max(0,q-ord))
         ERK_natural_ord(il,ib) = ord;
         ERK_natural_ordred(il,ib) = max(0,q-ord);

       end
     end
   end
   
   % solve with 'forced' approach
   if (run_forced)
     fprintf('\n  Results with forced approach:\n')

     % set problem-defining functions  (t must be scalar in each)
     fn = @f_timedep_bdry2;
     fe = @fe_timedep_bdry2;
     fi = @fi_timedep_bdry2;
     Jn = @J_timedep_bdry2;
     Ji = @J_timedep_bdry2;
     bdry = @enforce_timedep_bdry2;

     % run DIRK tests
     if (length(DIRKmethods) > 0)
       for ib = 1:length(DIRKmethods)

         mname = DIRKmethods{ib};
         B = butcher(mname);  s = numel(B(1,:))-1;  q = B(s+1,1);
         fprintf('\n    DIRK integrator: %s (order = %i)\n',mname,q)
         errs = zeros(size(hvals));
         order = zeros(length(hvals)-1,1);
         for ih = 1:length(hvals)
           hval = hvals(ih)/lambda;
           fprintf('      h = %.5e,',hval);
           [t,Y,ns,nl,cf,af] = solve_DIRK_bdry(fn, Jn, bdry, tout, Y0, B, rtol, atol, hval, hval);
           errs(ih) = sqrt(sum(sum((Y-Yref).^2))/numel(Y));
           if (isnan(errs(ih)) || isinf(errs(ih)) || errs(ih) > 10)
             errs(ih) = 1;
           end
           DIRK_forced_rmserrs(il,ib,ih) = errs(ih);
           if (ih>1)
             order(ih-1) = log( errs(ih)/errs(ih-1) ) / log( hvals(ih)/hvals(ih-1) );
             fprintf('  rmserr = %.2e,  lsolves = %i,  order = %.2e\n',errs(ih),nl,order(ih-1));
           else
             fprintf('  rmserr = %.2e,  lsolves = %i\n',errs(ih),nl);
           end
           fprintf('      newton conv fails = %i, temporal error fails = %i\n',cf,af);
         end
         s = sort(order);  ord = s(end-1);
         fprintf('    estimated order = %g,  reduction = %g\n',ord,max(0,q-ord))
         DIRK_forced_ord(il,ib) = ord;
         DIRK_forced_ordred(il,ib) = max(0,q-ord);

       end
     end

     % run ARK tests
     if ((length(ARKImethods) == length(ARKEmethods)) && (length(ARKImethods)>0))
       for ib = 1:length(ARKImethods)

         mname1 = ARKEmethods{ib};
         Be = butcher(mname1);  s = numel(Be(1,:))-1;  q = Be(s+1,1);
         mname2 = ARKImethods{ib};
         Bi = butcher(mname2);
         fprintf('\n    ARK integrator: %s/%s (order = %i)\n',mname1,mname2,q)
         errs = zeros(size(hvals));
         order = zeros(length(hvals)-1,1);
         for ih = 1:length(hvals)
           hval = hvals(ih)/lambda;
           fprintf('      h = %.5e,',hval);
           [t,Y,ns,nl,cf,af] = solve_ARK_bdry(fe, fi, Ji, bdry, tout, Y0, Be, Bi, rtol, atol, hval, hval);
           errs(ih) = sqrt(sum(sum((Y-Yref).^2))/numel(Y));
           if (isnan(errs(ih)) || isinf(errs(ih)) || errs(ih) > 10)
             errs(ih) = 1;
           end
           ARK_forced_rmserrs(il,ib,ih) = errs(ih);
           if (ih>1)
             order(ih-1) = log( errs(ih)/errs(ih-1) ) / log( hvals(ih)/hvals(ih-1) );
             fprintf('  rmserr = %.2e,  lsolves = %i,  order = %.2e\n',errs(ih),nl,order(ih-1));
           else
             fprintf('  rmserr = %.2e,  lsolves = %i\n',errs(ih),nl);
           end
           fprintf('      newton conv fails = %i, temporal error fails = %i\n',cf,af);
         end
         s = sort(order);  ord = s(end-1);
         fprintf('    estimated order = %g,  reduction = %g\n',ord,max(0,q-ord))
         ARK_forced_ord(il,ib) = ord;
         ARK_forced_ordred(il,ib) = max(0,q-ord);

       end
     end

     % run ERK tests
     if ((length(ERKmethods) > 0) && (abs(lambda)<=100))
       for ib = 1:length(ERKmethods)
         mname = ERKmethods{ib};
         B = butcher(mname);  s = numel(B(1,:))-1;  q = B(s+1,1);
         fprintf('\n   ERK integrator: %s (order = %i)\n',mname,q)
         errs_rms = zeros(size(hvals));
         errs_max = zeros(size(hvals));
         order = zeros(length(hvals)-1,1);
         for ih = 1:length(hvals)
           hval = hvals(ih)/lambda/100;
           fprintf('   h = %.5e,',hval);
           [t,Y,ns,af] = solve_ERK_bdry(fn, Es, bdry, tout, Y0, B, rtol, atol, hval, hval);
           errs(ih) = sqrt(sum(sum((Y-Yref).^2))/numel(Y));
           if (isnan(errs(ih)) || isinf(errs(ih)) || errs(ih) > 10)
             errs(ih) = 1;
           end
           ARK_forced_rmserrs(il,ib,ih) = errs(ih);
           if (ih>1)
             order(ih-1) = log( errs(ih)/errs(ih-1) ) / log( hvals(ih)/hvals(ih-1) );
             fprintf('  rmserr = %.2e,  lsolves = %i,  order = %.2e\n',errs(ih),nl,order(ih-1));
           else
             fprintf('  rmserr = %.2e,  lsolves = %i\n',errs(ih),nl);
           end
           fprintf('      temporal error fails = %i\n',af);
         end
         s = sort(order);  ord = s(end-1);
         fprintf('    estimated order = %g,  reduction = %g\n',ord,max(0,q-ord))
         ERK_forced_ord(il,ib) = ord;
         ERK_forced_ordred(il,ib) = max(0,q-ord);

       end
     end
   end
   

   % solve with 'approximate' approach
   if (run_approx)
     fprintf('\n  Results with approx approach:\n')

     % set problem-defining functions  (t must be scalar in each)
     fn = @f_timedep_bdry3;
     fe = @fe_timedep_bdry3;
     fi = @fi_timedep_bdry;
     Jn = @J_timedep_bdry;
     Ji = @J_timedep_bdry;
     global bDotApprox

     % run DIRK tests
     if (length(DIRKmethods) > 0)
       for ib = 1:length(DIRKmethods)

         mname = DIRKmethods{ib};
         B = butcher(mname);  s = numel(B(1,:))-1;  q = B(s+1,1);
         fprintf('\n    DIRK integrator: %s (order = %i)\n',mname,q)
         errs = zeros(size(hvals));
         order = zeros(length(hvals)-1,1);
         for ih = 1:length(hvals)
           hval = hvals(ih)/lambda;
           fprintf('      h = %.5e,',hval);
           bDotApprox.nmax = max(3,q);
           bDotApprox.nstored = 0;
           bDotApprox.b1 = zeros(bDotApprox.nmax,1);
           bDotApprox.b2 = zeros(bDotApprox.nmax,1);
           bDotApprox.t = zeros(bDotApprox.nmax,1);
           bDotApprox.h = hval;
           [t,Y,ns,nl,cf,af] = solve_DIRK(fn, Jn, tout, Y0, B, rtol, atol, hval, hval, 0);
           errs(ih) = sqrt(sum(sum((Y-Yref).^2))/numel(Y));
           if (isnan(errs(ih)) || isinf(errs(ih)) || errs(ih) > 10)
             errs(ih) = 1;
           end
           DIRK_approx_rmserrs(il,ib,ih) = errs(ih);
           if (ih>1)
             order(ih-1) = log( errs(ih)/errs(ih-1) ) / log( hvals(ih)/hvals(ih-1) );
             fprintf('  rmserr = %.2e,  lsolves = %i,  order = %.2e\n',errs(ih),nl,order(ih-1));
           else
             fprintf('  rmserr = %.2e,  lsolves = %i\n',errs(ih),nl);
           end
           fprintf('      newton conv fails = %i, temporal error fails = %i\n',cf,af);
         end
         s = sort(order);  ord = s(end-1);
         fprintf('    estimated order = %g,  reduction = %g\n',ord,max(0,q-ord))
         DIRK_approx_ord(il,ib) = ord;
         DIRK_approx_ordred(il,ib) = max(0,q-ord);

       end
     end

     % run ARK tests
     if ((length(ARKImethods) == length(ARKEmethods)) && (length(ARKImethods)>0))
       for ib = 1:length(ARKImethods)

         mname1 = ARKEmethods{ib};
         Be = butcher(mname1);  s = numel(Be(1,:))-1;  q = Be(s+1,1);
         mname2 = ARKImethods{ib};
         Bi = butcher(mname2);
         fprintf('\n    ARK integrator: %s/%s (order = %i)\n',mname1,mname2,q)
         errs = zeros(size(hvals));
         order = zeros(length(hvals)-1,1);
         for ih = 1:length(hvals)
           hval = hvals(ih)/lambda;
           fprintf('      h = %.5e,',hval);
           bDotApprox.nmax = max(3,q);
           bDotApprox.nstored = 0;
           bDotApprox.b1 = zeros(bDotApprox.nmax,1);
           bDotApprox.b2 = zeros(bDotApprox.nmax,1);
           bDotApprox.t = zeros(bDotApprox.nmax,1);
           bDotApprox.h = hval;
           [t,Y,ns,nl,cf,af] = solve_ARK(fe, fi, Ji, tout, Y0, Be, Bi, rtol, atol, hval, hval, 0);
           errs(ih) = sqrt(sum(sum((Y-Yref).^2))/numel(Y));
           if (isnan(errs(ih)) || isinf(errs(ih)) || errs(ih) > 10)
             errs(ih) = 1;
           end
           ARK_approx_rmserrs(il,ib,ih) = errs(ih);
           if (ih>1)
             order(ih-1) = log( errs(ih)/errs(ih-1) ) / log( hvals(ih)/hvals(ih-1) );
             fprintf('  rmserr = %.2e,  lsolves = %i,  order = %.2e\n',errs(ih),nl,order(ih-1));
           else
             fprintf('  rmserr = %.2e,  lsolves = %i\n',errs(ih),nl);
           end
           fprintf('      newton conv fails = %i, temporal error fails = %i\n',cf,af);
         end
         s = sort(order);  ord = s(end-1);
         fprintf('    estimated order = %g,  reduction = %g\n',ord,max(0,q-ord))
         ARK_approx_ord(il,ib) = ord;
         ARK_approx_ordred(il,ib) = max(0,q-ord);

       end
     end

     % run ERK tests
     if ((length(ERKmethods) > 0) && (abs(lambda)<=100))
       for ib = 1:length(ERKmethods)
         mname = ERKmethods{ib};
         B = butcher(mname);  s = numel(B(1,:))-1;  q = B(s+1,1);
         fprintf('\n   ERK integrator: %s (order = %i)\n',mname,q)
         errs_rms = zeros(size(hvals));
         errs_max = zeros(size(hvals));
         order = zeros(length(hvals)-1,1);
         for ih = 1:length(hvals)
           hval = hvals(ih)/lambda/20;
           fprintf('      h = %.5e,',hval);
           bDotApprox.nmax = max(3,q);
           bDotApprox.nstored = 0;
           bDotApprox.b1 = zeros(bDotApprox.nmax,1);
           bDotApprox.b2 = zeros(bDotApprox.nmax,1);
           bDotApprox.t = zeros(bDotApprox.nmax,1);
           bDotApprox.h = hval;
           [t,Y,ns,af] = solve_ERK(fn, Es, tout, Y0, B, rtol, atol, hval, hval);
           errs(ih) = sqrt(sum(sum((Y-Yref).^2))/numel(Y));
           if (isnan(errs(ih)) || isinf(errs(ih)) || errs(ih) > 10)
             errs(ih) = 1;
           end
           ARK_approx_rmserrs(il,ib,ih) = errs(ih);
           if (ih>1)
             order(ih-1) = log( errs(ih)/errs(ih-1) ) / log( hvals(ih)/hvals(ih-1) );
             fprintf('  rmserr = %.2e,  lsolves = %i,  order = %.2e\n',errs(ih),nl,order(ih-1));
           else
             fprintf('  rmserr = %.2e,  lsolves = %i\n',errs(ih),nl);
           end
           fprintf('      temporal error fails = %i\n',af);
         end
         s = sort(order);  ord = s(end-1);
         fprintf('    estimated order = %g,  reduction = %g\n',ord,max(0,q-ord))
         ERK_approx_ord(il,ib) = ord;
         ERK_approx_ordred(il,ib) = max(0,q-ord);

       end
     end
   end
   
end


% store statistics to disk
if (store_output)
  save timedep_bdry_data.mat;
end

% stop diary
if (create_diary)
  diary off
end

% end of script
