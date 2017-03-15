% driver for Van der Pol ODE test problem:
%    u' = v
%    v' = (v - v*u^2)/ep - u
% where u(0) = 2,  v(0) = 0, and ep = 0.2, integrated over 
% the time interval [0,12].
%
% Daniel R. Reynolds
% Department of Mathematics
% Southern Methodist University
% March 2017
% All Rights Reserved
clear

% set problem parameters
fn = @f_p1;
fe = @fe_p1;
fi = @fi_p1;
Jn = @J_p1;
Ji = @Ji_p1;
Es = @EStab_p1;
Tf = 12;
tout = linspace(0,Tf,100);
hmin = 1e-6;
hmax = 1.0;
rtol = 1e-3;
atol = 1e-14*ones(2,1);
global Pdata;
Pdata.ep = 0.2;
u0 = 2;
v0 = 0;
Y0 = [u0; v0];

% plot "true" solution
opts = odeset('RelTol',1e-12, 'AbsTol',atol,'InitialStep',hmin/10, 'MaxStep',hmax);
[t,Ytrue] = ode15s(fn, tout, Y0, opts);
figure()
plot(tout,Ytrue)
xlabel('t','FontSize',12), ylabel('y','FontSize',12)
title('Van der Pol test','FontSize',14)
set(gca,'FontSize',12)
print('-dpng','vanderPol')


% run with a diagonally-implicit RK method
mname = 'TRBDF2-ESDIRK';
B = butcher(mname);  s = numel(B(1,:))-1;
fprintf('\nRunning with DIRK integrator: %s (order = %i)\n',mname,B(s+1,1))
[t,Y,ns,nl] = solve_DIRK(fn, Jn, tout, Y0, B, rtol, atol, hmin, hmax);
err_max = max(max(abs(Y'-Ytrue)));
err_rms = sqrt(sum(sum((Y'-Ytrue).^2))/numel(Y));
fprintf('Accuracy/Work Results:\n')
fprintf('   maxerr = %.5e,  rmserr = %.5e\n',err_max, err_rms);
fprintf('   steps = %i (stages = %i), linear solves = %i\n',ns,ns*s,nl);


% run with an additive RK method
mname1 = 'ARK3(2)4L[2]SA-ERK';
Be = butcher(mname1);  s = numel(Be(1,:))-1;
mname2 = 'ARK3(2)4L[2]SA-ESDIRK';
Bi = butcher(mname2);
fprintf('\nRunning with ARK integrator: %s/%s (order = %i)\n',...
        mname1,mname2,Be(s+1,1))
[t,Y,ns,nl] = solve_ARK(fe, fi, Ji, tout, Y0, Be, Bi, rtol, atol, hmin, hmax);
err_max = max(max(abs(Y'-Ytrue)));
err_rms = sqrt(sum(sum((Y'-Ytrue).^2))/numel(Y));
fprintf('Accuracy/Work Results:\n')
fprintf('   maxerr = %.5e,  rmserr = %.5e\n',err_max, err_rms);
fprintf('   steps = %i (stages = %i), linear solves = %i\n',ns,ns*s,nl);


% run with an explicit RK method
mname = 'Bogacki-Shampine-ERK';
B = butcher(mname);  s = numel(B(1,:))-1;
fprintf('\nRunning with ERK integrator: %s (order = %i)\n',mname,B(s+1,1))
[t,Y,ns] = solve_ERK(fn, Es, tout, Y0, B, rtol, atol, hmin, hmax);
err_max = max(max(abs(Y'-Ytrue)));
err_rms = sqrt(sum(sum((Y'-Ytrue).^2))/numel(Y));
fprintf('Accuracy/Work Results:\n')
fprintf('   maxerr = %.5e,  rmserr = %.5e\n',err_max, err_rms);
fprintf('   steps = %i (stages = %i)\n',ns,ns*s);


% end of script
