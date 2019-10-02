% driver for a simple test problem with time-dependent
% mass "matrix" and analytical solution,
%    M(t) dy/dt = (1+t^2)*lamda*(y-atan(t)) + 1, -3<t<7
%    y(0) = atan(-3)
% where we use the simple mass matrix
%    M(t) = 1+t^2
% This has analytical solution
%    y(t) = atan(t).
% The stiffness of the problem is directly proportional to the
% value of "lamda"; since this test problem is used to ascertain
% order of accuracy (in the absence of order reduction), we use
% a relatively benign value of lamda=-10.
%
% This program solves the problem with either DIRK, ARK and ERK
% methods.  Unlike other test problems, this one uses a variety 
% of fixed step sizes to estimate the order of convergence for 
% the method (to determine whether the mass matrix time dependence 
% is handled appropriately).  If the estimated convergence rate is
% too far afield from expectations, we set "maxerr" and "rmserr"
% to large values so that regression tests are marked as "failed".
%
% Daniel R. Reynolds
% Department of Mathematics
% Southern Methodist University
% March 2017
% All Rights Reserved
clear

% set problem parameters
T0 = -3;
Tf = 7;
tout = linspace(T0,Tf,101);
hvals = [0.1, 0.05, 0.025, 0.01, 0.005, 0.0025, 0.001];
rtol = 1e-3;
atol = 1e-14*ones(2,1);
lambda = -10;
Y0 = atan(T0);


% set problem-defining functions
Mn = @(t)   1+t.^2;
fn = @(t,y) Mn(t)*lambda*(y - atan(t)) + 1;
fe = @(t,y) 1 - Mn(t)*lambda*atan(t);
fi = @(t,y) Mn(t)*lambda*y;
Jn = @(t,y) Mn(t)*lambda;
Ji = @(t,y) Mn(t)*lambda;
Es = @(t,y) 1/Mn(t)/abs(lambda);
ytrue = @(t) atan(t);


% perform tests with a diagonally-implicit RK method, algorithm 0
mname = 'Kvaerno(5,3,4)-ESDIRK';
B = butcher(mname);  s = numel(B(1,:))-1;
fprintf('\nRunning with DIRK integrator, algorithm 0: %s (order = %i)\n',mname,B(s+1,1))
errs_rms = zeros(size(hvals));
errs_max = zeros(size(hvals));
order = zeros(length(hvals)-1,1);
for ih = 1:length(hvals)
   fprintf('   h = %.5e,',hvals(ih));
   [t,Y,ns,nl] = solve_DIRK_mass(Mn, fn, Jn, tout, Y0, B, rtol, atol, hvals(ih), hvals(ih), 0);
   errs_max(ih) = max(max(abs(Y-ytrue(t))));
   errs_rms(ih) = sqrt(sum(sum((Y-ytrue(t)).^2))/numel(Y));
   if (ih>1) 
      order(ih-1) = log( errs_rms(ih)/errs_rms(ih-1) ) / log( hvals(ih)/hvals(ih-1) );
   end
   fprintf('   maxerr = %.5e,   rmserr = %.5e\n',errs_max(ih), errs_rms(ih));
   fprintf('   steps = %i (stages = %i), linear solves = %i\n',ns,ns*s,nl);
end
fprintf('Order of accuracy estimates (based on RMS errors above):\n')
for ih = 2:length(hvals)
   fprintf('  %g', order(ih-1) )
end
fprintf('\nOverall order of accuracy estimate = %g\n\n',sum(order)/length(order))


% perform tests with a diagonally-implicit RK method, algorithm 1
mname = 'Kvaerno(5,3,4)-ESDIRK';
B = butcher(mname);  s = numel(B(1,:))-1;
fprintf('\nRunning with DIRK integrator, algorithm 1: %s (order = %i)\n',mname,B(s+1,1))
errs_rms = zeros(size(hvals));
errs_max = zeros(size(hvals));
order = zeros(length(hvals)-1,1);
for ih = 1:length(hvals)
   fprintf('   h = %.5e,',hvals(ih));
   [t,Y,ns,nl] = solve_DIRK_mass(Mn, fn, Jn, tout, Y0, B, rtol, atol, hvals(ih), hvals(ih), 1);
   errs_max(ih) = max(max(abs(Y-ytrue(t))));
   errs_rms(ih) = sqrt(sum(sum((Y-ytrue(t)).^2))/numel(Y));
   if (ih>1) 
      order(ih-1) = log( errs_rms(ih)/errs_rms(ih-1) ) / log( hvals(ih)/hvals(ih-1) );
   end
   fprintf('   maxerr = %.5e,   rmserr = %.5e\n',errs_max(ih), errs_rms(ih));
   fprintf('   steps = %i (stages = %i), linear solves = %i\n',ns,ns*s,nl);
end
fprintf('Order of accuracy estimates (based on RMS errors above):\n')
for ih = 2:length(hvals)
   fprintf('  %g', order(ih-1) )
end
fprintf('\nOverall order of accuracy estimate = %g\n\n',sum(order)/length(order))


% run with an additive RK method, algorithm 0
mname1 = 'ARK4(3)7L[2]SA-ERK';
Be = butcher(mname1);  s = numel(Be(1,:))-1;
mname2 = 'ARK4(3)7L[2]SA-ESDIRK';
Bi = butcher(mname2);
fprintf('\nRunning with ARK integrator, algorithm 0: %s/%s (order = %i)\n',mname1,mname2,Be(s+1,1))
errs_rms = zeros(size(hvals));
errs_max = zeros(size(hvals));
order = zeros(length(hvals)-1,1);
for ih = 1:length(hvals)
   fprintf('   h = %.5e,',hvals(ih));
   [t,Y,ns,nl] = solve_ARK_mass(Mn, fe, fi, Ji, tout, Y0, Be, Bi, rtol, atol, hvals(ih), hvals(ih), 0);
   errs_max(ih) = max(max(abs(Y-ytrue(t))));
   errs_rms(ih) = sqrt(sum(sum((Y-ytrue(t)).^2))/numel(Y));
   if (ih>1) 
      order(ih-1) = log( errs_rms(ih)/errs_rms(ih-1) ) / log( hvals(ih)/hvals(ih-1) );
   end
   fprintf('   maxerr = %.5e,   rmserr = %.5e\n',errs_max(ih), errs_rms(ih));
   fprintf('   steps = %i (stages = %i), linear solves = %i\n',ns,ns*s,nl);
end
fprintf('Order of accuracy estimates (based on RMS errors above):\n')
for ih = 2:length(hvals)
   fprintf('  %g', order(ih-1) )
end
fprintf('\nOverall order of accuracy estimate = %g\n\n',sum(order)/length(order))


% run with an additive RK method, algorithm 1
mname1 = 'ARK4(3)7L[2]SA-ERK';
Be = butcher(mname1);  s = numel(Be(1,:))-1;
mname2 = 'ARK4(3)7L[2]SA-ESDIRK';
Bi = butcher(mname2);
fprintf('\nRunning with ARK integrator, algorithm 1: %s/%s (order = %i)\n',mname1,mname2,Be(s+1,1))
errs_rms = zeros(size(hvals));
errs_max = zeros(size(hvals));
order = zeros(length(hvals)-1,1);
for ih = 1:length(hvals)
   fprintf('   h = %.5e,',hvals(ih));
   [t,Y,ns,nl] = solve_ARK_mass(Mn, fe, fi, Ji, tout, Y0, Be, Bi, rtol, atol, hvals(ih), hvals(ih), 1);
   errs_max(ih) = max(max(abs(Y-ytrue(t))));
   errs_rms(ih) = sqrt(sum(sum((Y-ytrue(t)).^2))/numel(Y));
   if (ih>1) 
      order(ih-1) = log( errs_rms(ih)/errs_rms(ih-1) ) / log( hvals(ih)/hvals(ih-1) );
   end
   fprintf('   maxerr = %.5e,   rmserr = %.5e\n',errs_max(ih), errs_rms(ih));
   fprintf('   steps = %i (stages = %i), linear solves = %i\n',ns,ns*s,nl);
end
fprintf('Order of accuracy estimates (based on RMS errors above):\n')
for ih = 2:length(hvals)
   fprintf('  %g', order(ih-1) )
end
fprintf('\nOverall order of accuracy estimate = %g\n\n',sum(order)/length(order))


% run with an explicit RK method, algorithm 0
mname = 'ERK-4-4';
B = butcher(mname);  s = numel(B(1,:))-1;
fprintf('\nRunning with ERK integrator, algorithm 0: %s (order = %i)\n',mname,B(s+1,1))
errs_rms = zeros(size(hvals));
errs_max = zeros(size(hvals));
order = zeros(length(hvals)-1,1);
for ih = 1:length(hvals)
   fprintf('   h = %.5e,',hvals(ih));
   [t,Y,ns] = solve_ERK_mass(Mn, fn, Es, tout, Y0, B, rtol, atol, hvals(ih), hvals(ih), 0);
   errs_max(ih) = max(max(abs(Y-ytrue(t))));
   errs_rms(ih) = sqrt(sum(sum((Y-ytrue(t)).^2))/numel(Y));
   if (ih>1) 
      order(ih-1) = log( errs_rms(ih)/errs_rms(ih-1) ) / log( hvals(ih)/hvals(ih-1) );
   end
   fprintf('   maxerr = %.5e,   rmserr = %.5e\n',errs_max(ih), errs_rms(ih));
   fprintf('   steps = %i (stages = %i)\n',ns,ns*s);
end
fprintf('Order of accuracy estimates (based on RMS errors above):\n')
for ih = 2:length(hvals)
   fprintf('  %g', order(ih-1) )
end
fprintf('\nOverall order of accuracy estimate = %g\n\n',sum(order)/length(order))


% run with an explicit RK method, algorithm 1
mname = 'ERK-4-4';
B = butcher(mname);  s = numel(B(1,:))-1;
fprintf('\nRunning with ERK integrator, algorithm 1: %s (order = %i)\n',mname,B(s+1,1))
errs_rms = zeros(size(hvals));
errs_max = zeros(size(hvals));
order = zeros(length(hvals)-1,1);
for ih = 1:length(hvals)
   fprintf('   h = %.5e,',hvals(ih));
   [t,Y,ns] = solve_ERK_mass(Mn, fn, Es, tout, Y0, B, rtol, atol, hvals(ih), hvals(ih), 1);
   errs_max(ih) = max(max(abs(Y-ytrue(t))));
   errs_rms(ih) = sqrt(sum(sum((Y-ytrue(t)).^2))/numel(Y));
   if (ih>1) 
      order(ih-1) = log( errs_rms(ih)/errs_rms(ih-1) ) / log( hvals(ih)/hvals(ih-1) );
   end
   fprintf('   maxerr = %.5e,   rmserr = %.5e\n',errs_max(ih), errs_rms(ih));
   fprintf('   steps = %i (stages = %i)\n',ns,ns*s);
end
fprintf('Order of accuracy estimates (based on RMS errors above):\n')
for ih = 2:length(hvals)
   fprintf('  %g', order(ih-1) )
end
fprintf('\nOverall order of accuracy estimate = %g\n\n',sum(order)/length(order))


% end of script
