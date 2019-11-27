% Driver for heat equation with time-dependent boundary conditions:
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

% set test flags
plot_ref = false;
surf_vs_contour = 'surf';
true_error = false;

% set problem parameters
lambda = 1;
xl = 0;
xr = pi;
m = 50;
Tf = 5/lambda;
tout = linspace(0,Tf,11);
hvals = [0.5, 0.25, 0.1, 0.05, 0.025, 0.01, 0.005]/lambda;
rtol = 1e-7;
atol = 1e-13*ones(m,1);
global Pdata;
Pdata.c = [1, 2, 3, 4, 5];
Pdata.lambda = lambda;
Pdata.m = m;
Pdata.xspan = linspace(xl,xr,m)';
Pdata.dx = (Pdata.xspan(end)-Pdata.xspan(1))/(m-1);

% initial conditions
Y0 = zeros(size(Pdata.xspan));

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

% set true solution  (t must be scalar, x must be a column vector; hard-coded for 5 modes)
utrue = @(t,x) ( cos(x*1)*(Pdata.c(1))/(Pdata.lambda*1)*(1-exp(-Pdata.lambda*t)) ...
               + cos(x*2)*(Pdata.c(2))/(Pdata.lambda*4)*(1-exp(-Pdata.lambda*t*4)) ...
               + cos(x*3)*(Pdata.c(3))/(Pdata.lambda*9)*(1-exp(-Pdata.lambda*t*9)) ...
               + cos(x*4)*(Pdata.c(4))/(Pdata.lambda*16)*(1-exp(-Pdata.lambda*t*16)) ...
               + cos(x*5)*(Pdata.c(5))/(Pdata.lambda*25)*(1-exp(-Pdata.lambda*t*25)) );

% generate true solution
Ytrue = zeros(m,length(tout));
for i=1:length(tout)
   Ytrue(:,i) = utrue(tout(i),Pdata.xspan);
end

% general output header information
fprintf('\nTime-dependent boundary test problem:\n')
fprintf('    spatial domain: [%g, %g]\n',Pdata.xspan(1),Pdata.xspan(end))
fprintf('    time domain: [0, %g]\n',Tf)
fprintf('    spatial mesh nodes: %i\n',m)
if (true_error)
   fprintf('    using true solution for accuracy estimates\n')
else
   fprintf('    using reference solution for accuracy estimates\n')
end


% generate reference solution with ode15s
if (plot_ref || ~true_error)
   fprintf('\nRunning with ode15s:\n')
   opts = odeset('RelTol',1e-13, 'AbsTol',atol,'Jacobian',Jn);
   [t,Yref] = ode15s(fn, tout, Y0, opts);
   Yref = Yref';
   err_max = max(max(abs(Yref-Ytrue)));
   err_rms = sqrt(sum(sum((Yref-Ytrue).^2))/numel(Ytrue));
   fprintf('  Reference solution accuracy:')
   fprintf('  maxerr = %.5e,  rmserr = %.5e\n',err_max, err_rms);
end
if (plot_ref)
   if (strcmp(surf_vs_contour,'surf'))
      subplot(1,3,1), surf(Pdata.xspan,t,Yref'), shading interp
      subplot(1,3,2), surf(Pdata.xspan,t,Ytrue'), shading interp
      subplot(1,3,3), surf(Pdata.xspan,t,Yref'-Ytrue'), shading interp
   else
      subplot(1,3,1), contourf(Pdata.xspan,t,Yref'), colorbar
      subplot(1,3,2), contourf(Pdata.xspan,t,Ytrue'), colorbar
      subplot(1,3,3), contourf(Pdata.xspan,t,Yref'-Ytrue'), colorbar
   end
   subplot(1,3,1), title('ODE15s','FontSize',12), xlabel('x','FontSize',12), ylabel('t','FontSize',12)
   subplot(1,3,2), title('True','FontSize',12), xlabel('x','FontSize',12), ylabel('t','FontSize',12)
   subplot(1,3,3), title('Error','FontSize',12), xlabel('x','FontSize',12), ylabel('t','FontSize',12)
   set(gca,'FontSize',12)
   drawnow
end

% run with a diagonally-implicit RK method, algorithm 0
mname = 'Kvaerno(7,4,5)-ESDIRK';
B = butcher(mname);  s = numel(B(1,:))-1;
fprintf('\nRunning with DIRK integrator, algorithm 0: %s (order = %i)\n',mname,B(s+1,1))
errs_rms = zeros(size(hvals));
errs_max = zeros(size(hvals));
order = zeros(length(hvals)-1,1);
for ih = 1:length(hvals)
   h = hvals(ih);
   fprintf('   h = %.5e,',h);
   [t,Y,ns,nl,cf,af] = solve_DIRK(fn, Jn, tout, Y0, B, rtol, atol, h, h, 0);
   if (true_error)
      errs_max(ih) = max(max(abs(Y-Ytrue)));
      errs_rms(ih) = sqrt(sum(sum((Y-Ytrue).^2))/numel(Y));
   else
      errs_max(ih) = max(max(abs(Y-Yref)));
      errs_rms(ih) = sqrt(sum(sum((Y-Yref).^2))/numel(Y));
   end
   if (ih>1) 
      order(ih-1) = log( errs_rms(ih)/errs_rms(ih-1) ) / log( hvals(ih)/hvals(ih-1) );
      fprintf(' maxerr = %.5e,   rmserr = %.5e,   order = %.5e\n',errs_max(ih), errs_rms(ih), order(ih-1));
   else
      fprintf(' maxerr = %.5e,   rmserr = %.5e\n',errs_max(ih), errs_rms(ih));
   end
   fprintf('   steps = %i (stages = %i), linear solves = %i\n',ns,ns*s,nl);
   fprintf('   newton conv fails = %i, temporal error fails = %i\n',cf,af);
end
fprintf('Order of accuracy estimates (based on RMS errors above):\n')
for ih = 2:length(hvals)
   fprintf('  %g', order(ih-1) )
end
s=sort(order);
fprintf('\nOverall order of accuracy estimate = %g\n',sum(order)/length(order))
fprintf('Max 2 order of accuracy estimates = %g, %g\n\n',s(end-1:end))


% run with a diagonally-implicit RK method, algorithm 1
mname = 'Kvaerno(7,4,5)-ESDIRK';
B = butcher(mname);  s = numel(B(1,:))-1;
fprintf('\nRunning with DIRK integrator, algorithm 1: %s (order = %i)\n',mname,B(s+1,1))
errs_rms = zeros(size(hvals));
errs_max = zeros(size(hvals));
order = zeros(length(hvals)-1,1);
for ih = 1:length(hvals)
   h = hvals(ih);
   fprintf('   h = %.5e,',h);
   [t,Y,ns,nl,cf,af] = solve_DIRK(fn, Jn, tout, Y0, B, rtol, atol, h, h, 1);
   if (true_error)
      errs_max(ih) = max(max(abs(Y-Ytrue)));
      errs_rms(ih) = sqrt(sum(sum((Y-Ytrue).^2))/numel(Y));
   else
      errs_max(ih) = max(max(abs(Y-Yref)));
      errs_rms(ih) = sqrt(sum(sum((Y-Yref).^2))/numel(Y));
   end
   if (ih>1) 
      order(ih-1) = log( errs_rms(ih)/errs_rms(ih-1) ) / log( hvals(ih)/hvals(ih-1) );
      fprintf(' maxerr = %.5e,   rmserr = %.5e,   order = %.5e\n',errs_max(ih), errs_rms(ih), order(ih-1));
   else
      fprintf(' maxerr = %.5e,   rmserr = %.5e\n',errs_max(ih), errs_rms(ih));
   end
   fprintf('   steps = %i (stages = %i), linear solves = %i\n',ns,ns*s,nl);
   fprintf('   newton conv fails = %i, temporal error fails = %i\n',cf,af);
end
fprintf('Order of accuracy estimates (based on RMS errors above):\n')
for ih = 2:length(hvals)
   fprintf('  %g', order(ih-1) )
end
s=sort(order);
fprintf('\nOverall order of accuracy estimate = %g\n',sum(order)/length(order))
fprintf('Max 2 order of accuracy estimates = %g, %g\n\n',s(end-1:end))


% run with an ARK method, algorithm 0
mname1 = 'ARK5(4)8L[2]SA-ERK';
Be = butcher(mname1);  s = numel(Be(1,:))-1;
mname2 = 'ARK5(4)8L[2]SA-ESDIRK';
Bi = butcher(mname2);  s = numel(Bi(1,:))-1;
fprintf('\nRunning with ARK integrator, algorithm 0: %s/%s (order = %i)\n',mname1,mname2,Be(s+1,1))
errs_rms = zeros(size(hvals));
errs_max = zeros(size(hvals));
order = zeros(length(hvals)-1,1);
for ih = 1:length(hvals)
   h = hvals(ih);
   fprintf('   h = %.5e,',h);
   [t,Y,ns,nl,cf,af] = solve_ARK(fe, fi, Ji, tout, Y0, Be, Bi, rtol, atol, h, h, 0);
   if (true_error)
      errs_max(ih) = max(max(abs(Y-Ytrue)));
      errs_rms(ih) = sqrt(sum(sum((Y-Ytrue).^2))/numel(Y));
   else
      errs_max(ih) = max(max(abs(Y-Yref)));
      errs_rms(ih) = sqrt(sum(sum((Y-Yref).^2))/numel(Y));
   end
   if (ih>1) 
      order(ih-1) = log( errs_rms(ih)/errs_rms(ih-1) ) / log( hvals(ih)/hvals(ih-1) );
      fprintf(' maxerr = %.5e,   rmserr = %.5e,   order = %.5e\n',errs_max(ih), errs_rms(ih), order(ih-1));
   else
      fprintf(' maxerr = %.5e,   rmserr = %.5e\n',errs_max(ih), errs_rms(ih));
   end
   fprintf('   steps = %i (stages = %i), linear solves = %i\n',ns,ns*s,nl);
   fprintf('   newton conv fails = %i, temporal error fails = %i\n',cf,af);
end
fprintf('Order of accuracy estimates (based on RMS errors above):\n')
for ih = 2:length(hvals)
   fprintf('  %g', order(ih-1) )
end
s=sort(order);
fprintf('\nOverall order of accuracy estimate = %g\n',sum(order)/length(order))
fprintf('Max 2 order of accuracy estimates = %g, %g\n\n',s(end-1:end))


% run with an ARK method, algorithm 1
mname1 = 'ARK5(4)8L[2]SA-ERK';
Be = butcher(mname1);  s = numel(Be(1,:))-1;
mname2 = 'ARK5(4)8L[2]SA-ESDIRK';
Bi = butcher(mname2);  s = numel(Bi(1,:))-1;
fprintf('\nRunning with ARK integrator, algorithm 1: %s/%s (order = %i)\n',mname1,mname2,Be(s+1,1))
errs_rms = zeros(size(hvals));
errs_max = zeros(size(hvals));
order = zeros(length(hvals)-1,1);
for ih = 1:length(hvals)
   h = hvals(ih);
   fprintf('   h = %.5e,',h);
   [t,Y,ns,nl,cf,af] = solve_ARK(fe, fi, Ji, tout, Y0, Be, Bi, rtol, atol, h, h, 1);
   if (true_error)
      errs_max(ih) = max(max(abs(Y-Ytrue)));
      errs_rms(ih) = sqrt(sum(sum((Y-Ytrue).^2))/numel(Y));
   else
      errs_max(ih) = max(max(abs(Y-Yref)));
      errs_rms(ih) = sqrt(sum(sum((Y-Yref).^2))/numel(Y));
   end
   if (ih>1) 
      order(ih-1) = log( errs_rms(ih)/errs_rms(ih-1) ) / log( hvals(ih)/hvals(ih-1) );
      fprintf(' maxerr = %.5e,   rmserr = %.5e,   order = %.5e\n',errs_max(ih), errs_rms(ih), order(ih-1));
   else
      fprintf(' maxerr = %.5e,   rmserr = %.5e\n',errs_max(ih), errs_rms(ih));
   end
   fprintf('   steps = %i (stages = %i), linear solves = %i\n',ns,ns*s,nl);
   fprintf('   newton conv fails = %i, temporal error fails = %i\n',cf,af);
end
fprintf('Order of accuracy estimates (based on RMS errors above):\n')
for ih = 2:length(hvals)
   fprintf('  %g', order(ih-1) )
end
s=sort(order);
fprintf('\nOverall order of accuracy estimate = %g\n',sum(order)/length(order))
fprintf('Max 2 order of accuracy estimates = %g, %g\n\n',s(end-1:end))


% run with an explicit RK method
mname = 'Fehlberg-ERK';
B = butcher(mname);  s = numel(B(1,:))-1;
fprintf('\nRunning with ERK integrator: %s (order = %i)\n',mname,B(s+1,1))
errs_rms = zeros(size(hvals));
errs_max = zeros(size(hvals));
order = zeros(length(hvals)-1,1);
for ih = 1:length(hvals)
   h = hvals(ih)/100;
   fprintf('   h = %.5e,',h);
   [t,Y,ns,af] = solve_ERK(fn, Es, tout, Y0, B, rtol, atol, h, h);
   if (true_error)
      errs_max(ih) = max(max(abs(Y-Ytrue)));
      errs_rms(ih) = sqrt(sum(sum((Y-Ytrue).^2))/numel(Y));
   else
      errs_max(ih) = max(max(abs(Y-Yref)));
      errs_rms(ih) = sqrt(sum(sum((Y-Yref).^2))/numel(Y));
   end
   if (ih>1) 
      order(ih-1) = log( errs_rms(ih)/errs_rms(ih-1) ) / log( hvals(ih)/hvals(ih-1) );
      fprintf(' maxerr = %.5e,   rmserr = %.5e,   order = %.5e\n',errs_max(ih), errs_rms(ih), order(ih-1));
   else
      fprintf(' maxerr = %.5e,   rmserr = %.5e\n',errs_max(ih), errs_rms(ih));
   end
   fprintf('   steps = %i (stages = %i)\n',ns,ns*s);
   fprintf('   temporal error fails = %i\n',af);
end
fprintf('Order of accuracy estimates (based on RMS errors above):\n')
for ih = 2:length(hvals)
   fprintf('  %g', order(ih-1) )
end
s=sort(order);
fprintf('\nOverall order of accuracy estimate = %g\n',sum(order)/length(order))
fprintf('Max 2 order of accuracy estimates = %g, %g\n\n',s(end-1:end))


% end of script
