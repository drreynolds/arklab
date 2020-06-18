% driver for stiff brusselator test problem:
%      u' = a - (w+1)*u + u^2*v,
%      v' = w*u - u^2*v,
%      w' = (b-w)/ep - u*w,
% where u(0) = 1.2, v(0) = 3.1 and w(0) = 3, with prameters a=1,
% b=3.5 and ep=5e-6.  We evaluate over the time interval [0,10].  
%
% Daniel R. Reynolds
% Department of Mathematics
% Southern Methodist University
% March 2017
% All Rights Reserved
clear

% set problem parameters
Tf = 10;
tout = linspace(0,Tf,101);
hmin = 1e-6;
hmax = 1.0;
rtol = 1e-6;
atol = 1e-14*ones(3,1);
a = 1; 
b = 3.5; 
ep = 1e-3;
u0 = 1.2;
v0 = 3.1;
w0 = 3;
Y0 = [u0; v0; w0];

% methods to test
test_sdirk  = false;
test_ark    = true;
test_erk    = false;
test_alg0   = true;
test_alg1   = false;

% generate plot of 'true' solution
make_plot   = false;

% set problem-defining functions
fn = @(t,y) [a - (y(3)+1)*y(1) + y(1)*y(1)*y(2); y(3)*y(1) - y(1)*y(1)*y(2); (b-y(3))/ep - y(3)*y(1)];
fe = @(t,y) [a - (y(3)+1)*y(1) + y(1)*y(1)*y(2); y(3)*y(1) - y(1)*y(1)*y(2); -y(3)*y(1)];
fi = @(t,y) [0; 0; (b-y(3))/ep];
Jn = @(t,y) [-(y(3)+1) + 2*y(1)*y(2), y(1)*y(1), -y(1);  y(3) - 2*y(1)*y(2), -y(1)*y(1), y(1); -y(3), 0, -1/ep - y(1)];
Ji = @(t,y) [0, 0, 0; 0, 0, 0; 0, 0, -1/ep];
Es = @(t,y) 1/max(abs(eig(Jn(t,y))));

% generate reference solution
opts = odeset('RelTol',1e-12, 'AbsTol',atol,'InitialStep',hmin/10, 'MaxStep',hmax);
[t,Ytrue] = ode15s(fn, tout, Y0, opts);
save('bruss_times.txt', '-ascii', 't')
save('bruss_solution.txt', '-ascii', 'Ytrue')

% plot reference solution
if (make_plot)
  figure()
  plot(tout,Ytrue)
  xlabel('t','FontSize',12), ylabel('y','FontSize',12)
  title('Brusselator ODE test','FontSize',14)
  set(gca,'FontSize',12)
  print('-dpng','brusselator')
end

% run with a diagonally-implicit RK method, algorithm 0
if (test_sdirk && test_alg0)
  mname = 'Cash(5,3,4)-SDIRK';
  B = butcher(mname);  s = numel(B(1,:))-1;
  fprintf('\nRunning with DIRK integrator, algorithm 0: %s (order = %i)\n',mname,B(s+1,1))
  [t,Y,ns,nl,cf,af] = solve_DIRK(fn, Jn, tout, Y0, B, rtol, atol, hmin, hmax, 0);
  err_max = max(max(abs(Y'-Ytrue)));
  err_rms = sqrt(sum(sum((Y'-Ytrue).^2))/numel(Y));
  fprintf('Accuracy/Work Results:\n')
  fprintf('   maxerr = %.5e,  rmserr = %.5e\n',err_max, err_rms);
  fprintf('   steps = %i (stages = %i), linear solves = %i\n',ns,ns*s,nl);
  fprintf('   newton conv fails = %i, temporal error fails = %i\n',cf,af);
end

% run with a diagonally-implicit RK method, algorithm 1
if (test_sdirk && test_alg1)
  mname = 'Cash(5,3,4)-SDIRK';
  B = butcher(mname);  s = numel(B(1,:))-1;
  fprintf('\nRunning with DIRK integrator, algorithm 1: %s (order = %i)\n',mname,B(s+1,1))
  [t,Y,ns,nl,cf,af] = solve_DIRK(fn, Jn, tout, Y0, B, rtol, atol, hmin, hmax, 1);
  err_max = max(max(abs(Y'-Ytrue)));
  err_rms = sqrt(sum(sum((Y'-Ytrue).^2))/numel(Y));
  fprintf('Accuracy/Work Results:\n')
  fprintf('   maxerr = %.5e,  rmserr = %.5e\n',err_max, err_rms);
  fprintf('   steps = %i (stages = %i), linear solves = %i\n',ns,ns*s,nl);
  fprintf('   newton conv fails = %i, temporal error fails = %i\n',cf,af);
end

% run with an ARK method, algorithm 0
if (test_ark && test_alg0)
  %mname1 = 'ARK4(3)6L[2]SA-ERK';
  mname1 = 'ARK3(2)4L[2]SA-ERK';
  Be = butcher(mname1);  s = numel(Be(1,:))-1;
  %mname2 = 'ARK4(3)6L[2]SA-ESDIRK';
  mname2 = 'ARK3(2)4L[2]SA-ESDIRK';
  Bi = butcher(mname2);
  fprintf('\nRunning with ARK integrator, algorithm 0: %s/%s (order = %i)\n',mname1,mname2,Be(s+1,1))
  keyboard
  [t,Y,ns,nl,cf,af] = solve_ARK(fe, fi, Ji, tout, Y0, Be, Bi, rtol, atol, hmin, hmax, 0);
  err_max = max(max(abs(Y'-Ytrue)));
  err_rms = sqrt(sum(sum((Y'-Ytrue).^2))/numel(Y));
  fprintf('Accuracy/Work Results:\n')
  fprintf('   maxerr = %.5e,  rmserr = %.5e\n',err_max, err_rms);
  fprintf('   steps = %i (stages = %i), linear solves = %i\n',ns,ns*s,nl);
  fprintf('   newton conv fails = %i, temporal error fails = %i\n',cf,af);
end

% run with an ARK method, algorithm 1
if (test_ark && test_alg1)
  %mname1 = 'ARK4(3)6L[2]SA-ERK';
  mname1 = 'ARK3(2)4L[2]SA-ERK';
  Be = butcher(mname1);  s = numel(Be(1,:))-1;
  %mname2 = 'ARK4(3)6L[2]SA-ESDIRK';
  mname2 = 'ARK3(2)4L[2]SA-ESDIRK';
  Bi = butcher(mname2);
  fprintf('\nRunning with ARK integrator, algorithm 1: %s/%s (order = %i)\n',mname1,mname2,Be(s+1,1))
  [t,Y,ns,nl,cf,af] = solve_ARK(fe, fi, Ji, tout, Y0, Be, Bi, rtol, atol, hmin, hmax, 1);
  err_max = max(max(abs(Y'-Ytrue)));
  err_rms = sqrt(sum(sum((Y'-Ytrue).^2))/numel(Y));
  fprintf('Accuracy/Work Results:\n')
  fprintf('   maxerr = %.5e,  rmserr = %.5e\n',err_max, err_rms);
  fprintf('   steps = %i (stages = %i), linear solves = %i\n',ns,ns*s,nl);
  fprintf('   newton conv fails = %i, temporal error fails = %i\n',cf,af);
end

% run with an explicit RK method
if (test_erk)
  mname = 'Fehlberg-ERK';
  B = butcher(mname);  s = numel(B(1,:))-1;
  fprintf('\nRunning with ERK integrator: %s (order = %i)\n',mname,B(s+1,1))
  [t,Y,ns,af] = solve_ERK(fn, Es, tout, Y0, B, rtol, atol, hmin, hmax);
  err_max = max(max(abs(Y'-Ytrue)));
  err_rms = sqrt(sum(sum((Y'-Ytrue).^2))/numel(Y));
  fprintf('Accuracy/Work Results:\n')
  fprintf('   maxerr = %.5e,  rmserr = %.5e\n',err_max, err_rms);
  fprintf('   steps = %i (stages = %i)\n',ns,ns*s);
  fprintf('   temporal error fails = %i\n',af);
end

% end of script
