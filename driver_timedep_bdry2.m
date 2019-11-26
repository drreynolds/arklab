% Driver for heat equation with time-dependent Dirichlet boundary conditions:
%      u_t = lambda*u_xx + f,  0<x<pi/2,  0<t<4/lambda
%      u(0,x) = 0,             0<x<pi/2,
%      u(t,0) = b1(t),                    0<t<4/lambda
%      u(t,pi) = b2(t),                   0<t<4/lambda
% The forcing and boundary terms have the form
%    f(t,x) = \sum_{k=1}^M c_k \cos(kx)
%    b1(t)  = \sum_{k=1}^M c_k/(lambda*k^2) (1 - e^{-t*lambda*k^2})
%    b2(t)  = -\sum_{k=1}^M c_k/(lambda*k^2) (1 - e^{-t*lambda*k^2})
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

% set problem parameters
lambda = 1;
xl = 0;
xr = pi;
m = 50;
Tf = 4/lambda;
tout = linspace(0,Tf,10);
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
fn = @f_timedep_bdry2;
fe = @fe_timedep_bdry2;
fi = @fi_timedep_bdry2;
Jn = @J_timedep_bdry2;
Ji = @J_timedep_bdry2;
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
bdry = @enforce_timedep_bdry2;

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
fprintf('\nTime-dependent boundary test problem 2:\n')
fprintf('    spatial domain: [%g, %g]\n',Pdata.xspan(1),Pdata.xspan(end))
fprintf('    time domain: [0, %g]\n',Tf)
fprintf('    spatial mesh nodes: %i\n',m)


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
   [t,Y,ns,nl,cf,af] = solve_DIRK_bdry(fn, Jn, bdry, tout, Y0, B, rtol, atol, h, h, 0);
   errs_max(ih) = max(max(abs(Y-Ytrue)));
   errs_rms(ih) = sqrt(sum(sum((Y-Ytrue).^2))/numel(Y));
   if (ih>1) 
      order(ih-1) = log( errs_rms(ih)/errs_rms(ih-1) ) / log( hvals(ih)/hvals(ih-1) );
   end
   fprintf('   maxerr = %.5e,   rmserr = %.5e\n',errs_max(ih), errs_rms(ih));
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


% $$$ % run with a diagonally-implicit RK method, algorithm 1
% $$$ mname = 'Kvaerno(7,4,5)-ESDIRK';
% $$$ B = butcher(mname);  s = numel(B(1,:))-1;
% $$$ fprintf('\nRunning with DIRK integrator, algorithm 1: %s (order = %i)\n',mname,B(s+1,1))
% $$$ errs_rms = zeros(size(hvals));
% $$$ errs_max = zeros(size(hvals));
% $$$ order = zeros(length(hvals)-1,1);
% $$$ for ih = 1:length(hvals)
% $$$    h = hvals(ih);
% $$$    fprintf('   h = %.5e,',h);
% $$$    [t,Y,ns,nl,cf,af] = solve_DIRK_bdry(fn, Jn, bdry, tout, Y0, B, rtol, atol, h, h, 1);
% $$$    errs_max(ih) = max(max(abs(Y-Ytrue)));
% $$$    errs_rms(ih) = sqrt(sum(sum((Y-Ytrue).^2))/numel(Y));
% $$$    if (ih>1) 
% $$$       order(ih-1) = log( errs_rms(ih)/errs_rms(ih-1) ) / log( hvals(ih)/hvals(ih-1) );
% $$$    end
% $$$    fprintf('   maxerr = %.5e,   rmserr = %.5e\n',errs_max(ih), errs_rms(ih));
% $$$    fprintf('   steps = %i (stages = %i), linear solves = %i\n',ns,ns*s,nl);
% $$$    fprintf('   newton conv fails = %i, temporal error fails = %i\n',cf,af);
% $$$ end
% $$$ fprintf('Order of accuracy estimates (based on RMS errors above):\n')
% $$$ for ih = 2:length(hvals)
% $$$    fprintf('  %g', order(ih-1) )
% $$$ end
% $$$ s=sort(order);
% $$$ fprintf('\nOverall order of accuracy estimate = %g\n',sum(order)/length(order))
% $$$ fprintf('Max 2 order of accuracy estimates = %g, %g\n\n',s(end-1:end))


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
   [t,Y,ns,nl,cf,af] = solve_ARK_bdry(fe, fi, Ji, bdry, tout, Y0, Be, Bi, rtol, atol, h, h, 0);
   errs_max(ih) = max(max(abs(Y-Ytrue)));
   errs_rms(ih) = sqrt(sum(sum((Y-Ytrue).^2))/numel(Y));
   if (ih>1) 
      order(ih-1) = log( errs_rms(ih)/errs_rms(ih-1) ) / log( hvals(ih)/hvals(ih-1) );
   end
   fprintf('   maxerr = %.5e,   rmserr = %.5e\n',errs_max(ih), errs_rms(ih));
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


% $$$ % run with an ARK method, algorithm 1
% $$$ mname1 = 'ARK5(4)8L[2]SA-ERK';
% $$$ Be = butcher(mname1);  s = numel(Be(1,:))-1;
% $$$ mname2 = 'ARK5(4)8L[2]SA-ESDIRK';
% $$$ Bi = butcher(mname2);  s = numel(Bi(1,:))-1;
% $$$ fprintf('\nRunning with ARK integrator, algorithm 1: %s/%s (order = %i)\n',mname1,mname2,Be(s+1,1))
% $$$ errs_rms = zeros(size(hvals));
% $$$ errs_max = zeros(size(hvals));
% $$$ order = zeros(length(hvals)-1,1);
% $$$ for ih = 1:length(hvals)
% $$$    h = hvals(ih);
% $$$    fprintf('   h = %.5e,',h);
% $$$    [t,Y,ns,nl,cf,af] = solve_ARK_bdry(fe, fi, Ji, bdry, tout, Y0, Be, Bi, rtol, atol, h, h, 1);
% $$$    errs_max(ih) = max(max(abs(Y-Ytrue)));
% $$$    errs_rms(ih) = sqrt(sum(sum((Y-Ytrue).^2))/numel(Y));
% $$$    if (ih>1) 
% $$$       order(ih-1) = log( errs_rms(ih)/errs_rms(ih-1) ) / log( hvals(ih)/hvals(ih-1) );
% $$$    end
% $$$    fprintf('   maxerr = %.5e,   rmserr = %.5e\n',errs_max(ih), errs_rms(ih));
% $$$    fprintf('   steps = %i (stages = %i), linear solves = %i\n',ns,ns*s,nl);
% $$$    fprintf('   newton conv fails = %i, temporal error fails = %i\n',cf,af);
% $$$ end
% $$$ fprintf('Order of accuracy estimates (based on RMS errors above):\n')
% $$$ for ih = 2:length(hvals)
% $$$    fprintf('  %g', order(ih-1) )
% $$$ end
% $$$ s=sort(order);
% $$$ fprintf('\nOverall order of accuracy estimate = %g\n',sum(order)/length(order))
% $$$ fprintf('Max 2 order of accuracy estimates = %g, %g\n\n',s(end-1:end))


% run with an explicit RK method
mname = 'Merson-5-4-ERK';
B = butcher(mname);  s = numel(B(1,:))-1;
fprintf('\nRunning with ERK integrator: %s (order = %i)\n',mname,B(s+1,1))
errs_rms = zeros(size(hvals));
errs_max = zeros(size(hvals));
order = zeros(length(hvals)-1,1);
for ih = 1:length(hvals)
   h = hvals(ih)/100;
   fprintf('   h = %.5e,',h);
   [t,Y,ns,af] = solve_ERK_bdry(fn, Es, bdry, tout, Y0, B, rtol, atol, h, h);
   errs_max(ih) = max(max(abs(Y-Ytrue)));
   errs_rms(ih) = sqrt(sum(sum((Y-Ytrue).^2))/numel(Y));
   if (ih>1) 
      order(ih-1) = log( errs_rms(ih)/errs_rms(ih-1) ) / log( hvals(ih)/hvals(ih-1) );
   end
   fprintf('   maxerr = %.5e,   rmserr = %.5e\n',errs_max(ih), errs_rms(ih));
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
