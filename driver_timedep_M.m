% Driver for a simple test problem with time-dependent
% mass "matrix" and analytical solution,
%    M(t)*dy/dt = M(t)*lamda*(y-g(t)) + M(t)*g'(t), t0<=t<=tf
%    y(0) = g(t0)
% where we use the simple mass "matrix"
%    M(t) = gamma/g'(t)
% under the condition that g'(t) ~= 0 for all t0<=t<=tf.
% This has analytical solution
%    y(t) = g(t).
% Here we use g(t) = atan(t) and time interval [-3,7],
% and hence
%    g'(t) = 1/(1+t^2),
%    M(t) = gamma*(1+t^2)
% The stiffness of the problem is directly proportional to the
% value of "lamda".  The 'units' of the mass matrix are
% proportional to the value of "gamma".
%
% Note: to get a 'baseline' value of order reduction for each RK
% method, we also run each method on the above problem using M(t)=1.
%
% This program solves the problem with either DIRK, ARK and ERK
% methods.  Unlike other test problems, this one uses a variety 
% of fixed step sizes to estimate the order of convergence for 
% the method (to determine whether the mass matrix time dependence 
% is handled appropriately).
%
% Daniel R. Reynolds
% Department of Mathematics
% Southern Methodist University
% October 2019
% All Rights Reserved
clear

% start output diary
!\rm timedep_M_results.txt
diary timedep_M_results.txt

% set problem parameters
lambdas = [-10, -100, -1000, -1e4, -1e5];
gammas = [1, 10, 100];
g  = @(t) atan(t);
gp = @(t) 1/(1+t.^2);
T0 = -3;
Tf = 7;

% set testing parameters
tout = linspace(T0,Tf,101);
hvals = [0.05, 0.025, 0.01, 0.005, 0.0025, 0.001, 0.0005, 0.00025];
rtol = 1e-5;     % since problem is linear in y, and tests use
atol = 1e-14;    % fixed stepsizes, tolerances are essentially
                 % irrelevant since 1 Newton iteration is 'exact'
Y0 = g(T0);
algs = [0, 1];
DIRKmethods = {'SDIRK-2-2','EDIRK-3-3','Kvaerno(5,3,4)-ESDIRK','Kvaerno(7,4,5)-ESDIRK'};
ARKEmethods = {'ARK(2,3,2)-ERK','Ascher(2,3,3)-ERK','ARK4(3)7L[2]SA-ERK','ARK5(4)8L[2]SA-ERK'};
ARKImethods = {'ARK(2,3,2)-SDIRK','Ascher(2,3,3)-SDIRK','ARK4(3)7L[2]SA-ESDIRK','ARK5(4)8L[2]SA-ESDIRK'};
ERKmethods  = {'ERK-1-1','ERK-2-2','ARK(2,3,2)-ERK','ERK-4-4','Dormand-Prince-ERK'};


% output general testing information
fprintf('\nRunning time-dependent mass matrix tests:\n');
fprintf('   Solver tolerances:  rtol = %g,  atol = %g\n', rtol, atol);


% create 'statistics' storage
DIRK_rmserrs = zeros(length(gammas),length(lambdas),length(DIRKmethods),length(algs),length(hvals));
DIRK_lsolves = zeros(length(gammas),length(lambdas),length(DIRKmethods),length(algs),length(hvals));
DIRK_ord = zeros(length(gammas),length(lambdas),length(DIRKmethods),length(algs));
DIRK_ordred = zeros(length(gammas),length(lambdas),length(DIRKmethods),length(algs));

ARK_rmserrs = zeros(length(gammas),length(lambdas),length(ARKEmethods),length(algs),length(hvals));
ARK_lsolves = zeros(length(gammas),length(lambdas),length(ARKEmethods),length(algs),length(hvals));
ARK_ord = zeros(length(gammas),length(lambdas),length(ARKEmethods),length(algs));
ARK_ordred = zeros(length(gammas),length(lambdas),length(ARKEmethods),length(algs));

ERK_rmserrs = zeros(length(gammas),length(lambdas),length(ERKmethods),length(algs),length(hvals));
ERK_ord = zeros(length(gammas),length(lambdas),length(ERKmethods),length(algs));
ERK_ordred = zeros(length(gammas),length(lambdas),length(ERKmethods),length(algs));


% loop over stiffness and mass matrix parameters
for ig = 1:length(gammas)
   gamma = gammas(ig);
   for il = 1:length(lambdas)
      lambda = lambdas(il);

      % set problem-defining functions
      Mn = @(t)   gamma/gp(t);
      fn = @(t,y) Mn(t)*lambda*(y - g(t)) + Mn(t)*gp(t);
      fe = @(t,y) Mn(t)*gp(t) - Mn(t)*lambda*g(t);
      fi = @(t,y) Mn(t)*lambda*y;
      Jn = @(t,y) Mn(t)*lambda;
      Ji = @(t,y) Mn(t)*lambda;
      Es = @(t,y) 1/Mn(t)/abs(lambda);

      % output general testing information
      fprintf('\n Stiffness factor lambda  = %g\n', lambda);
      fprintf(' Mass matrix factor gamma = %g\n', gamma);

      
      % run DIRK tests
      if (length(DIRKmethods) > 0) 
         for ib = 1:length(DIRKmethods)
            for ia = 1:length(algs)
               alg = algs(ia);
               mname = DIRKmethods{ib};
               B = butcher(mname);  s = numel(B(1,:))-1;  q = B(s+1,1);
               fprintf('\n  DIRK integrator, algorithm %i: %s (order = %i)\n',alg,mname,q)
               errs = zeros(size(hvals));
               order = zeros(length(hvals)-1,1);
               for ih = 1:length(hvals)
                  h = hvals(ih);
                  fprintf('    h = %.5e,',h);
                  [t,Y,ns,nl] = solve_DIRK_mass(Mn, fn, Jn, tout, Y0, B, rtol, atol, h, h, alg);
                  errs(ih) = sqrt(sum(sum((Y-g(t)).^2))/numel(Y));
                  DIRK_rmserrs(ig,il,ib,ia,ih) = errs(ih);
                  DIRK_lsolves(ig,il,ib,ia,ih) = nl;
                  if (ih>1) 
                     order(ih-1) = log( errs(ih)/errs(ih-1) ) / log( hvals(ih)/hvals(ih-1) );
                  end
                  fprintf('  rmserr = %.2e,  lsolves = %i\n',errs(ih),nl);
               end
               s = sort(order);  ord = s(end-1);
               fprintf('  estimated order = %g,  reduction = %g\n',ord,max(0,q-ord))
               DIRK_ord(ig,il,ib,ia) = ord;
               DIRK_ordred(ig,il,ib,ia,ih) = max(0,q-ord);
            end
         end
      end

      
      % run ARK tests
      if ((length(ARKImethods) == length(ARKEmethods)) && (length(ARKImethods)>0))
         for ib = 1:length(ARKImethods)
            for ia = 1:length(algs)
               alg = algs(ia);
               mname1 = ARKEmethods{ib};
               Be = butcher(mname1);  s = numel(Be(1,:))-1;  q = Be(s+1,1);
               mname2 = ARKImethods{ib};
               Bi = butcher(mname2);
               fprintf('\n  ARK integrator, algorithm %i: %s/%s (order = %i)\n',alg,mname1,mname2,q)
               errs = zeros(size(hvals));
               order = zeros(length(hvals)-1,1);
               for ih = 1:length(hvals)
                  h = hvals(ih);
                  fprintf('     h = %.5e,',h);
                  [t,Y,ns,nl] = solve_ARK_mass(Mn, fe, fi, Ji, tout, Y0, Be, Bi, rtol, atol, h, h, alg);
                  errs(ih) = sqrt(sum(sum((Y-g(t)).^2))/numel(Y));
                  ARK_rmserrs(ig,il,ib,ia,ih) = errs(ih);
                  ARK_lsolves(ig,il,ib,ia,ih) = nl;
                  if (ih>1) 
                     order(ih-1) = log( errs(ih)/errs(ih-1) ) / log( hvals(ih)/hvals(ih-1) );
                  end
                  fprintf('  rmserr = %.2e,  lsolves = %i\n',errs(ih),nl);
               end
               s = sort(order);  ord = s(end-1);
               fprintf('  estimated order = %g,  reduction = %g\n',ord,max(0,q-ord))
               ARK_ord(ig,il,ib,ia) = ord;
               ARK_ordred(ig,il,ib,ia,ih) = max(0,q-ord);
            end
         end
      end


      % run ERK tests
      if ((length(ERKmethods) > 0) && (abs(lambda)<=100))
         for ib = 1:length(ERKmethods)
            for ia = 1:length(algs)
               alg = algs(ia);
               mname = ERKmethods{ib};
               B = butcher(mname);  s = numel(B(1,:))-1;  q = B(s+1,1);
               fprintf('\n  ERK integrator, algorithm %i: %s (order = %i)\n',alg,mname,q)
               errs = zeros(size(hvals));
               order = zeros(length(hvals)-1,1);
               for ih = 1:length(hvals)
                  h = hvals(ih)/abs(lambda);
                  fprintf('     h = %.5e,',h);
                  [t,Y,ns] = solve_ERK_mass(Mn, fn, Es, tout, Y0, B, rtol, atol, h, h, alg);
                  errs(ih) = sqrt(sum(sum((Y-g(t)).^2))/numel(Y));
                  ERK_rmserrs(ig,il,ib,ia,ih) = errs(ih);
                  if (ih>1) 
                     order(ih-1) = log( errs(ih)/errs(ih-1) ) / log( hvals(ih)/hvals(ih-1) );
                  end
                  fprintf('  rmserr = %.2e\n',errs(ih));
               end
               s = sort(order);  ord = s(end-1);
               fprintf('  estimated order = %g,  reduction = %g\n',ord,max(0,q-ord))
               ERK_ord(ig,il,ib,ia) = ord;
               ERK_ordred(ig,il,ib,ia,ih) = max(0,q-ord);
            end
         end
      end

   
   end
end




% output general testing information
fprintf('\nRe-running tests without mass matrix (for baseline order reduction):\n');
fprintf('   Solver tolerances:  rtol = %g,  atol = %g\n', rtol, atol);


% create 'statistics' storage
DIRK_rmserrs_noM = zeros(length(lambdas),length(DIRKmethods),length(algs),length(hvals));
DIRK_lsolves_noM = zeros(length(lambdas),length(DIRKmethods),length(algs),length(hvals));
DIRK_ord_noM = zeros(length(lambdas),length(DIRKmethods),length(algs));
DIRK_ordred_noM = zeros(length(lambdas),length(DIRKmethods),length(algs));

ARK_rmserrs_noM = zeros(length(lambdas),length(ARKEmethods),length(algs),length(hvals));
ARK_lsolves_noM = zeros(length(lambdas),length(ARKEmethods),length(algs),length(hvals));
ARK_ord_noM = zeros(length(lambdas),length(ARKEmethods),length(algs));
ARK_ordred_noM = zeros(length(lambdas),length(ARKEmethods),length(algs));

ERK_rmserrs_noM = zeros(length(lambdas),length(ERKmethods),length(algs),length(hvals));
ERK_ord_noM = zeros(length(lambdas),length(ERKmethods),length(algs));
ERK_ordred_noM = zeros(length(lambdas),length(ERKmethods),length(algs));


% loop over stiffness parameters
for il = 1:length(lambdas)
   lambda = lambdas(il);

   % set problem-defining functions
   Mn = @(t)   1;
   fn = @(t,y) Mn(t)*lambda*(y - g(t)) + Mn(t)*gp(t);
   fe = @(t,y) Mn(t)*gp(t) - Mn(t)*lambda*g(t);
   fi = @(t,y) Mn(t)*lambda*y;
   Jn = @(t,y) Mn(t)*lambda;
   Ji = @(t,y) Mn(t)*lambda;
   Es = @(t,y) 1/Mn(t)/abs(lambda);
   
   % output general testing information
   fprintf('\n Stiffness factor lambda  = %g\n', lambda);
  
   
   % run DIRK tests
   if (length(DIRKmethods) > 0) 
      for ib = 1:length(DIRKmethods)
         for ia = 1:length(algs)
            alg = algs(ia);
            mname = DIRKmethods{ib};
            B = butcher(mname);  s = numel(B(1,:))-1;  q = B(s+1,1);
            fprintf('\n  DIRK integrator, algorithm %i: %s (order = %i)\n',alg,mname,q)
            errs = zeros(size(hvals));
            order = zeros(length(hvals)-1,1);
            for ih = 1:length(hvals)
               h = hvals(ih);
               fprintf('    h = %.5e,',h);
               [t,Y,ns,nl] = solve_DIRK_mass(Mn, fn, Jn, tout, Y0, B, rtol, atol, h, h, alg);
               errs(ih) = sqrt(sum(sum((Y-g(t)).^2))/numel(Y));
               DIRK_rmserrs_noM(il,ib,ia,ih) = errs(ih);
               DIRK_lsolves_noM(il,ib,ia,ih) = nl;
               if (ih>1) 
                  order(ih-1) = log( errs(ih)/errs(ih-1) ) / log( hvals(ih)/hvals(ih-1) );
               end
               fprintf('  rmserr = %.2e,  lsolves = %i\n',errs(ih),nl);
            end
            s = sort(order);  ord = s(end-1);
            fprintf('  estimated order = %g,  reduction = %g\n',ord,max(0,q-ord))
            DIRK_ord_noM(il,ib,ia) = ord;
            DIRK_ordred_noM(il,ib,ia,ih) = max(0,q-ord);
         end
      end
   end

   
   % run ARK tests
   if ((length(ARKImethods) == length(ARKEmethods)) && (length(ARKImethods)>0))
      for ib = 1:length(ARKImethods)
         for ia = 1:length(algs)
            alg = algs(ia);
            mname1 = ARKEmethods{ib};
            Be = butcher(mname1);  s = numel(Be(1,:))-1;  q = Be(s+1,1);
            mname2 = ARKImethods{ib};
            Bi = butcher(mname2);
            fprintf('\n  ARK integrator, algorithm %i: %s/%s (order = %i)\n',alg,mname1,mname2,q)
            errs = zeros(size(hvals));
            order = zeros(length(hvals)-1,1);
            for ih = 1:length(hvals)
               h = hvals(ih);
               fprintf('     h = %.5e,',h);
               [t,Y,ns,nl] = solve_ARK_mass(Mn, fe, fi, Ji, tout, Y0, Be, Bi, rtol, atol, h, h, alg);
               errs(ih) = sqrt(sum(sum((Y-g(t)).^2))/numel(Y));
               ARK_rmserrs_noM(il,ib,ia,ih) = errs(ih);
               ARK_lsolves_noM(il,ib,ia,ih) = nl;
               if (ih>1) 
                  order(ih-1) = log( errs(ih)/errs(ih-1) ) / log( hvals(ih)/hvals(ih-1) );
               end
               fprintf('  rmserr = %.2e,  lsolves = %i\n',errs(ih),nl);
            end
            s = sort(order);  ord = s(end-1);
            fprintf('  estimated order = %g,  reduction = %g\n',ord,max(0,q-ord))
            ARK_ord_noM(il,ib,ia) = ord;
            ARK_ordred_noM(il,ib,ia,ih) = max(0,q-ord);
         end
      end
   end


   % run ERK tests
   if ((length(ERKmethods) > 0) && (abs(lambda)<=100))
      for ib = 1:length(ERKmethods)
         for ia = 1:length(algs)
            alg = algs(ia);
            mname = ERKmethods{ib};
            B = butcher(mname);  s = numel(B(1,:))-1;  q = B(s+1,1);
            fprintf('\n  ERK integrator, algorithm %i: %s (order = %i)\n',alg,mname,q)
            errs = zeros(size(hvals));
            order = zeros(length(hvals)-1,1);
            for ih = 1:length(hvals)
               h = hvals(ih)/abs(lambda);
               fprintf('     h = %.5e,',h);
               [t,Y,ns] = solve_ERK_mass(Mn, fn, Es, tout, Y0, B, rtol, atol, h, h, alg);
               errs(ih) = sqrt(sum(sum((Y-g(t)).^2))/numel(Y));
               ERK_rmserrs_noM(il,ib,ia,ih) = errs(ih);
               if (ih>1) 
                  order(ih-1) = log( errs(ih)/errs(ih-1) ) / log( hvals(ih)/hvals(ih-1) );
               end
               fprintf('  rmserr = %.2e\n',errs(ih));
            end
            s = sort(order);  ord = s(end-1);
            fprintf('  estimated order = %g,  reduction = %g\n',ord,max(0,q-ord))
            ERK_ord_noM(il,ib,ia) = ord;
            ERK_ordred_noM(il,ib,ia,ih) = max(0,q-ord);
         end
      end
   end

   
end

% store statistics to disk
save timedep_M_data.mat;

% stop diary
diary off

% end of script
