% driver for brusselator PDE test problem:
%      u' = d1*u_xx + a - (b+1)*u + u^2*v,
%      v' = d2*v_xx + b*u - u^2*v,
% with parameters a=0.6, b=2, d1=d2=0.25, over the spatial 
% interval [0,1] and the time interval [0,80].  Initial 
% conditions are u(x,0) = a + 0.1*sin(pi*x), 
% v(x,0) = b/a + 0.1*sin(pi*x), and boundary conditions are 
% homogeneous Dirichlet.  We use a mesh with 100 spatial zones.
%
% Daniel R. Reynolds
% Department of Mathematics
% Southern Methodist University
% March 2017
% All Rights Reserved
clear

% set problem parameters
m = 100;
Tf = 10;
tout = linspace(0,Tf,100);
hmin = 1e-6;
hmax = 1.0;
rtol = 1e-3;
atol = 1e-14*ones(2*m,1);
global Pdata;
Pdata.a = 0.6;
Pdata.b = 2;
Pdata.d1 = 0.025;
Pdata.d2 = 0.025;
Pdata.m = m;
Pdata.dx = 1/(m-1);

% initial conditions
xspan = linspace(0,1,m)';
u0 = Pdata.a + 0.1*sin(pi*xspan);
v0 = Pdata.b/Pdata.a + 0.1*sin(pi*xspan);
Y0 = [u0; v0];

% set problem-defining functions
fn = @f_p3;
fe = @fe_p3;
fi = @fi_p3;
Jn = @J_p3;
Ji = @Ji_p3;
Es = @(t,y) 1/max(abs(eigs(Jn(t,y))));

% plot "true" solution at end of run
opts = odeset('RelTol',1e-12, 'AbsTol',atol,'InitialStep',hmin/10, 'MaxStep',hmax);
[t,Ytrue] = ode15s(fn, tout, Y0, opts);
figure()
plot(xspan,Ytrue(end,1:m),xspan,Ytrue(end,m+1:2*m))
xlabel('t','FontSize',12), ylabel('y','FontSize',12)
title('1D Brusselator PDE test','FontSize',14)
set(gca,'FontSize',12)
print('-dpng','brusselator1D')


% run with a diagonally-implicit RK method, algorithm 0
mname = 'Kvaerno(7,4,5)-ESDIRK';
B = butcher(mname);  s = numel(B(1,:))-1;
fprintf('\nRunning with DIRK integrator, algorithm 0: %s (order = %i)\n',mname,B(s+1,1))
[t,Y,ns,nl] = solve_DIRK(fn, Jn, tout, Y0, B, rtol, atol, hmin, hmax, 0);
err_max = max(max(abs(Y'-Ytrue)));
err_rms = sqrt(sum(sum((Y'-Ytrue).^2))/numel(Y));
fprintf('Accuracy/Work Results:\n')
fprintf('   maxerr = %.5e,  rmserr = %.5e\n',err_max, err_rms);
fprintf('   steps = %i (stages = %i), linear solves = %i\n',ns,ns*s,nl);


% run with a diagonally-implicit RK method, algorithm 1
mname = 'Kvaerno(7,4,5)-ESDIRK';
B = butcher(mname);  s = numel(B(1,:))-1;
fprintf('\nRunning with DIRK integrator, algorithm 1: %s (order = %i)\n',mname,B(s+1,1))
[t,Y,ns,nl] = solve_DIRK(fn, Jn, tout, Y0, B, rtol, atol, hmin, hmax, 1);
err_max = max(max(abs(Y'-Ytrue)));
err_rms = sqrt(sum(sum((Y'-Ytrue).^2))/numel(Y));
fprintf('Accuracy/Work Results:\n')
fprintf('   maxerr = %.5e,  rmserr = %.5e\n',err_max, err_rms);
fprintf('   steps = %i (stages = %i), linear solves = %i\n',ns,ns*s,nl);


% run with an ARK method, algorithm 0
mname1 = 'ARK5(4)8L[2]SA-ERK';
Be = butcher(mname1);  s = numel(Be(1,:))-1;
mname2 = 'ARK5(4)8L[2]SA-ESDIRK';
Bi = butcher(mname2);  s = numel(Bi(1,:))-1;
fprintf('\nRunning with ARK integrator, algorithm 0: %s/%s (order = %i)\n',mname1,mname2,Be(s+1,1))
[t,Y,ns,nl] = solve_ARK(fe, fi, Ji, tout, Y0, Be, Bi, rtol, atol, hmin, hmax, 0);
err_max = max(max(abs(Y'-Ytrue)));
err_rms = sqrt(sum(sum((Y'-Ytrue).^2))/numel(Y));
fprintf('Accuracy/Work Results:\n')
fprintf('   maxerr = %.5e,  rmserr = %.5e\n',err_max, err_rms);
fprintf('   steps = %i (stages = %i), linear solves = %i\n',ns,ns*s,nl);


% run with an ARK method, algorithm 1
mname1 = 'ARK5(4)8L[2]SA-ERK';
Be = butcher(mname1);  s = numel(Be(1,:))-1;
mname2 = 'ARK5(4)8L[2]SA-ESDIRK';
Bi = butcher(mname2);  s = numel(Bi(1,:))-1;
fprintf('\nRunning with ARK integrator, algorithm 1: %s/%s (order = %i)\n',mname1,mname2,Be(s+1,1))
[t,Y,ns,nl] = solve_ARK(fe, fi, Ji, tout, Y0, Be, Bi, rtol, atol, hmin, hmax, 1);
err_max = max(max(abs(Y'-Ytrue)));
err_rms = sqrt(sum(sum((Y'-Ytrue).^2))/numel(Y));
fprintf('Accuracy/Work Results:\n')
fprintf('   maxerr = %.5e,  rmserr = %.5e\n',err_max, err_rms);
fprintf('   steps = %i (stages = %i), linear solves = %i\n',ns,ns*s,nl);


% $$$ % run with an explicit RK method
% $$$ mname = 'Merson-5-4-ERK';
% $$$ B = butcher(mname);  s = numel(B(1,:))-1;
% $$$ fprintf('\nRunning with ERK integrator: %s (order = %i)\n',mname,B(s+1,1))
% $$$ [t,Y,ns] = solve_ERK(fn, Es, tout, Y0, B, rtol, atol, hmin, hmax);
% $$$ err_max = max(max(abs(Y'-Ytrue)));
% $$$ err_rms = sqrt(sum(sum((Y'-Ytrue).^2))/numel(Y));
% $$$ fprintf('Accuracy/Work Results:\n')
% $$$ fprintf('   maxerr = %.5e,  rmserr = %.5e\n',err_max, err_rms);
% $$$ fprintf('   steps = %i (stages = %i)\n',ns,ns*s);


% end of script
