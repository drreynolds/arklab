% driver for stiff brusselator test problem:
%      y' = f(y),
% where y = [u, v, w]', and f = [fu, fv, fw]' with
%      fu = a - (w+1)*u + u^2*v,
%      fv = w*u - u^2*v,
%      fw = (b-w)/ep - u*w,
% where u(0) = 1.2, v(0) = 3.1 and w(0) = 3, with parameters a=1,
% b=3.5 and ep=5e-6.  We evaluate over the time interval [0,10].
%
% This problem is augmented with a fixed mass matrix,
%    M = 2*rand(3,3);
% that modifies the problem to become
%    M*y' = M*f(y)
%
% Daniel R. Reynolds
% Department of Mathematics
% Southern Methodist University
% March 2017
% All Rights Reserved
clear

% set problem parameters
Tf = 10;
tout = linspace(0,Tf,100);
hmin = 1e-7;
hmax = 1.0;
rtol = 1e-3;
atol = 1e-14*ones(3,1);
a = 1; 
b = 3.5; 
ep = 1e-3;
u0 = 1.2;
v0 = 3.1;
w0 = 3;
Y0 = [u0; v0; w0];

% set problem-defining functions
rng('default');   % reset random number generator
rng(1);
M = 20*rand(3,3);
fref = @(t,y) [a - (y(3)+1)*y(1) + y(1)*y(1)*y(2); y(3)*y(1) - y(1)*y(1)*y(2); (b-y(3))/ep - y(3)*y(1)];
fn = @(t,y) M*[a - (y(3)+1)*y(1) + y(1)*y(1)*y(2); y(3)*y(1) - y(1)*y(1)*y(2); (b-y(3))/ep - y(3)*y(1)];
fe = @(t,y) M*[a - (y(3)+1)*y(1) + y(1)*y(1)*y(2); y(3)*y(1) - y(1)*y(1)*y(2); -y(3)*y(1)];
fi = @(t,y) M*[0; 0; (b-y(3))/ep];
Jn = @(t,y) M*[-(y(3)+1) + 2*y(1)*y(2), y(1)*y(1), -y(1);  y(3) - 2*y(1)*y(2), -y(1)*y(1), y(1); -y(3), 0, -1/ep - y(1)];
Ji = @(t,y) M*[0, 0, 0; 0, 0, 0; 0, 0, -1/ep];
%Es = @(t,y) 1/max(abs(eig(M\Jn(t,y))));
Es = @(t,y) Tf;

% plot "true" solution 
opts = odeset('RelTol',1e-12, 'AbsTol',atol,'InitialStep',hmin/10, 'MaxStep',hmax);
[t,Ytrue] = ode15s(fref, tout, Y0, opts);
figure()
plot(tout,Ytrue)
xlabel('t','FontSize',12), ylabel('y','FontSize',12)
title('Brusselator ODE test','FontSize',14)
set(gca,'FontSize',12)
print('-dpng','brusselator')


% run with a diagonally-implicit RK method, algorithm 0
mname = 'Cash(5,3,4)-SDIRK';
B = butcher(mname);  s = numel(B(1,:))-1;
fprintf('\nRunning with DIRK integrator, algorithm 0: %s (order = %i)\n',mname,B(s+1,1))
[t,Y,ns,nl] = solve_DIRK_mass(M, fn, Jn, tout, Y0, B, rtol, atol, hmin, hmax, 0);
err_max = max(max(abs(Y'-Ytrue)));
err_rms = sqrt(sum(sum((Y'-Ytrue).^2))/numel(Y));
fprintf('Accuracy/Work Results:\n')
fprintf('   maxerr = %.5e,  rmserr = %.5e\n',err_max, err_rms);
fprintf('   steps = %i (stages = %i), linear solves = %i\n',ns,ns*s,nl);


% run with a diagonally-implicit RK method, algorithm 1
mname = 'Cash(5,3,4)-SDIRK';
B = butcher(mname);  s = numel(B(1,:))-1;
fprintf('\nRunning with DIRK integrator, algorithm 1: %s (order = %i)\n',mname,B(s+1,1))
[t,Y,ns,nl] = solve_DIRK_mass(M, fn, Jn, tout, Y0, B, rtol, atol, hmin, hmax, 1);
err_max = max(max(abs(Y'-Ytrue)));
err_rms = sqrt(sum(sum((Y'-Ytrue).^2))/numel(Y));
fprintf('Accuracy/Work Results:\n')
fprintf('   maxerr = %.5e,  rmserr = %.5e\n',err_max, err_rms);
fprintf('   steps = %i (stages = %i), linear solves = %i\n',ns,ns*s,nl);


% run with an ARK method, algorithm 0
mname1 = 'ARK4(3)6L[2]SA-ERK';
Be = butcher(mname1);  s = numel(Be(1,:))-1;
mname2 = 'ARK4(3)6L[2]SA-ESDIRK';
Bi = butcher(mname2);
fprintf('\nRunning with ARK integrator, algorithm 0: %s/%s (order = %i)\n',mname1,mname2,Be(s+1,1))
[t,Y,ns,nl] = solve_ARK_mass(M, fe, fi, Ji, tout, Y0, Be, Bi, rtol, atol, hmin, hmax, 0);
err_max = max(max(abs(Y'-Ytrue)));
err_rms = sqrt(sum(sum((Y'-Ytrue).^2))/numel(Y));
fprintf('Accuracy/Work Results:\n')
fprintf('   maxerr = %.5e,  rmserr = %.5e\n',err_max, err_rms);
fprintf('   steps = %i (stages = %i), linear solves = %i\n',ns,ns*s,nl);


% run with an ARK method, algorithm 1
mname1 = 'ARK4(3)6L[2]SA-ERK';
Be = butcher(mname1);  s = numel(Be(1,:))-1;
mname2 = 'ARK4(3)6L[2]SA-ESDIRK';
Bi = butcher(mname2);
fprintf('\nRunning with ARK integrator, algorithm 1: %s/%s (order = %i)\n',mname1,mname2,Be(s+1,1))
[t,Y,ns,nl] = solve_ARK_mass(M, fe, fi, Ji, tout, Y0, Be, Bi, rtol, atol, hmin, hmax, 1);
err_max = max(max(abs(Y'-Ytrue)));
err_rms = sqrt(sum(sum((Y'-Ytrue).^2))/numel(Y));
fprintf('Accuracy/Work Results:\n')
fprintf('   maxerr = %.5e,  rmserr = %.5e\n',err_max, err_rms);
fprintf('   steps = %i (stages = %i), linear solves = %i\n',ns,ns*s,nl);


% run with an explicit RK method, algorithm 0
mname = 'Fehlberg-ERK';
B = butcher(mname);  s = numel(B(1,:))-1;
fprintf('\nRunning with ERK integrator, algorithm 0: %s (order = %i)\n',mname,B(s+1,1))
[t,Y,ns] = solve_ERK_mass(M, fn, Es, tout, Y0, B, rtol, atol, hmin, hmax, 0);
err_max = max(max(abs(Y'-Ytrue)));
err_rms = sqrt(sum(sum((Y'-Ytrue).^2))/numel(Y));
fprintf('Accuracy/Work Results:\n')
fprintf('   maxerr = %.5e,  rmserr = %.5e\n',err_max, err_rms);
fprintf('   steps = %i (stages = %i)\n',ns,ns*s);


% run with an explicit RK method, algorithm 1
mname = 'Fehlberg-ERK';
B = butcher(mname);  s = numel(B(1,:))-1;
fprintf('\nRunning with ERK integrator, algorithm 1: %s (order = %i)\n',mname,B(s+1,1))
[t,Y,ns] = solve_ERK_mass(M, fn, Es, tout, Y0, B, rtol, atol, hmin, hmax, 1);
err_max = max(max(abs(Y'-Ytrue)));
err_rms = sqrt(sum(sum((Y'-Ytrue).^2))/numel(Y));
fprintf('Accuracy/Work Results:\n')
fprintf('   maxerr = %.5e,  rmserr = %.5e\n',err_max, err_rms);
fprintf('   steps = %i (stages = %i)\n',ns,ns*s);


% end of script
