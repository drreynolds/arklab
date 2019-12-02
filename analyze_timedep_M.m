% Driver to analyze output from the script 'driver_timedep_M.m'.
% The primary purpose of that script, and these analyses, is to
% determine whether there is any appreciable performance difference
% between approaches that solve for stages, z, or implicit RHS
% vectors, k.  We thus examine a variety of metrics to assess the
% performance of both methods.
%
% The secondary purpose is to determine whether we notice any
% performance difference of these solvers on problems where M=I and
% where M=M(t).
%
% Daniel R. Reynolds
% Department of Mathematics
% Southern Methodist University
% October 2019
% All Rights Reserved
clear

% start output diary
!\rm timedep_M_analysis.txt
diary timedep_M_analysis.txt

% load results:
%   Test setup:
%     lambdas = stiffness parameters
%     gammas = mass matrix units
%     hvals = time step sizes
%     algs = solve for z (0) vs solve for k (1)
%     DIRKmethods = DIRK Butcher tables (strings)
%     ARKEmethods = ARK explicit Butcher tables (strings)
%     ARKImethods = ARK implicit Butcher tables (strings)
%     ERKmethods = ERK Butcher tables (strings)
%   Results for M=M(t):
%     DIRK_rmserrs(gammas,lambdas,DIRKmethods,algs,hvals)
%     DIRK_lsolves(gammas,lambdas,DIRKmethods,algs,hvals)
%     DIRK_ord(gammas,lambdas,DIRKmethods,algs)
%     DIRK_ordred(gammas,lambdas,DIRKmethods,algs)
%     ARK_rmserrs(gammas,lambdas,DIRKmethods,algs,hvals)
%     ARK_lsolves(gammas,lambdas,DIRKmethods,algs,hvals)
%     ARK_ord(gammas,lambdas,DIRKmethods,algs)
%     ARK_ordred(gammas,lambdas,DIRKmethods,algs)
%     ERK_rmserrs(gammas,lambdas,DIRKmethods,algs,hvals)
%     ERK_ord(gammas,lambdas,DIRKmethods,algs)
%     ERK_ordred(gammas,lambdas,DIRKmethods,algs)
%   Results for M=I:
%     DIRK_rmserrs_noM(lambdas,DIRKmethods,algs,hvals)
%     DIRK_lsolves_noM(lambdas,DIRKmethods,algs,hvals)
%     DIRK_ord_noM(lambdas,DIRKmethods,algs)
%     DIRK_ordred_noM(lambdas,DIRKmethods,algs)
%     ARK_rmserrs_noM(lambdas,DIRKmethods,algs,hvals)
%     ARK_lsolves_noM(lambdas,DIRKmethods,algs,hvals)
%     ARK_ord_noM(lambdas,DIRKmethods,algs)
%     ARK_ordred_noM(lambdas,DIRKmethods,algs)
%     ERK_rmserrs_noM(lambdas,DIRKmethods,algs,hvals)
%     ERK_ord_noM(lambdas,DIRKmethods,algs)
%     ERK_ordred_noM(lambdas,DIRKmethods,algs)
load timedep_M_data.mat

% set shortcut variables for testing array lengths
nl = length(lambdas);
ng = length(gammas);
nh = length(hvals);
nDIRK = length(DIRKmethods);
nARK = length(ARKEmethods);
nERK = length(ERKmethods);


%--- alg 0 vs alg 1 ---%
fprintf('\n\nComparing solving for z against solving for k:\n');


% Solution order
fprintf('\nOrder of accuracy:\n');

for ib = 1:length(DIRKmethods)
   mname = DIRKmethods{ib};
   B = butcher(mname);  s = numel(B(1,:))-1;  q = B(s+1,1);
   fprintf('\n  DIRK integrator: %s (order = %i)\n',mname,q)
   z_data = DIRK_ord(:,:,ib,1);
   k_data = DIRK_ord(:,:,ib,2);
   fprintf('      (z): mean = %g, min = %g, max = %g, std = %g\n', ...
           mean(z_data,'all','omitnan'), ...
           min(z_data,[],'all','omitnan'), ...
           max(z_data,[],'all','omitnan'), ...
           std(z_data,0,'all','omitnan'))
   fprintf('      (k): mean = %g, min = %g, max = %g, std = %g\n', ...
           mean(k_data,'all','omitnan'), ...
           min(k_data,[],'all','omitnan'), ...
           max(k_data,[],'all','omitnan'), ...
           std(k_data,0,'all','omitnan'))
   
   figure()
   z = squeeze(mean(z_data,1));
   k = squeeze(mean(k_data,1));
   semilogx(abs(lambdas),z,abs(lambdas),k)
   legend('(z)','(k)')
   xlabel('Stiffness: |\lambda|'), ylabel('order')
   title(sprintf('Measured order vs stiffness (method = %s)',mname))
   saveas(gcf, sprintf('timedep_M-order_vs_stiffness-DIRK_b%i.png',ib))
   
   figure()
   z = squeeze(mean(z_data,2));
   k = squeeze(mean(k_data,2));
   semilogx(gammas,z,gammas,k)
   legend('(z)','(k)')
   xlabel('Mass units: \gamma'), ylabel('order')
   title(sprintf('Measured order vs mass units (method = %s)',mname))
   saveas(gcf, sprintf('timedep_M-order_vs_massunits-DIRK_b%i.png',ib))
end

for ib = 1:length(ARKImethods)
   mname1 = ARKEmethods{ib};
   Be = butcher(mname1);  s = numel(Be(1,:))-1;  q = Be(s+1,1);
   mname2 = ARKImethods{ib};
   Bi = butcher(mname2);
   fprintf('\n  ARK integrator: %s/%s (order = %i)\n',mname1,mname2,q)
   z_data = ARK_ord(:,:,ib,1);
   k_data = ARK_ord(:,:,ib,2);
   fprintf('      (z): mean = %g, min = %g, max = %g, std = %g\n', ...
           mean(z_data,'all','omitnan'), ...
           min(z_data,[],'all','omitnan'), ...
           max(z_data,[],'all','omitnan'), ...
           std(z_data,0,'all','omitnan'))
   fprintf('      (k): mean = %g, min = %g, max = %g, std = %g\n', ...
           mean(k_data,'all','omitnan'), ...
           min(k_data,[],'all','omitnan'), ...
           max(k_data,[],'all','omitnan'), ...
           std(k_data,0,'all','omitnan'))
   
   figure()
   z = squeeze(mean(z_data,1));
   k = squeeze(mean(k_data,1));
   semilogx(abs(lambdas),z,abs(lambdas),k)
   legend('(z)','(k)')
   xlabel('Stiffness: |\lambda|'), ylabel('order')
   title(sprintf('Measured order vs stiffness (method = %s)',mname1))
   saveas(gcf, sprintf('timedep_M-order_vs_stiffness-ARK_b%i.png',ib))
   
   figure()
   z = squeeze(mean(z_data,2));
   k = squeeze(mean(k_data,2));
   semilogx(gammas,z,gammas,k)
   legend('(z)','(k)')
   xlabel('Mass units: \gamma'), ylabel('order')
   title(sprintf('Measured order vs mass units (method = %s)',mname1))
   saveas(gcf, sprintf('timedep_M-order_vs_massunits-ARK_b%i.png',ib))
end

nlam = sum(abs(lambdas) <= 100);
for ib = 1:length(ERKmethods)
   mname = ERKmethods{ib};
   B = butcher(mname);  s = numel(B(1,:))-1;  q = B(s+1,1);
   fprintf('\n  ERK integrator: %s (order = %i)\n',mname,q)
   z_data = ERK_ord(:,1:nlam,ib,1);
   k_data = ERK_ord(:,1:nlam,ib,2);
   fprintf('      (z): mean = %g, min = %g, max = %g, std = %g\n', ...
           mean(z_data,'all','omitnan'), ...
           min(z_data,[],'all','omitnan'), ...
           max(z_data,[],'all','omitnan'), ...
           std(z_data,0,'all','omitnan'))
   fprintf('      (k): mean = %g, min = %g, max = %g, std = %g\n', ...
           mean(k_data,'all','omitnan'), ...
           min(k_data,[],'all','omitnan'), ...
           max(k_data,[],'all','omitnan'), ...
           std(k_data,0,'all','omitnan'))
   
   figure()
   z = squeeze(mean(z_data,1));
   k = squeeze(mean(k_data,1));
   semilogx(abs(lambdas(1:nlam)),z,abs(lambdas(1:nlam)),k)
   legend('(z)','(k)')
   xlabel('Stiffness: |\lambda|'), ylabel('order')
   title(sprintf('Measured order vs stiffness (method = %s)',mname))
   saveas(gcf, sprintf('timedep_M-order_vs_stiffness-ERK_b%i.png',ib))
   
   figure()
   z = squeeze(mean(z_data,2));
   k = squeeze(mean(k_data,2));
   semilogx(gammas,z,gammas,k)
   legend('(z)','(k)')
   xlabel('Mass units: \gamma'), ylabel('order')
   title(sprintf('Measured order vs mass units (method = %s)',mname))
   saveas(gcf, sprintf('timedep_M-order_vs_massunits-ERK_b%i.png',ib))
end



% Solution accuracy
fprintf('\nAccuracy at smallest step size:\n');

for ib = 1:length(DIRKmethods)
   mname = DIRKmethods{ib};
   B = butcher(mname);  s = numel(B(1,:))-1;  q = B(s+1,1);
   fprintf('\n  DIRK integrator: %s (order = %i)\n',mname,q)
   z_data = DIRK_rmserrs(:,:,ib,1,end);
   k_data = DIRK_rmserrs(:,:,ib,2,end);
   fprintf('      (z): mean = %g, min = %g, max = %g, std = %g\n', ...
           mean(z_data,'all','omitnan'), ...
           min(z_data,[],'all','omitnan'), ...
           max(z_data,[],'all','omitnan'), ...
           std(z_data,0,'all','omitnan'))
   fprintf('      (k): mean = %g, min = %g, max = %g, std = %g\n', ...
           mean(k_data,'all','omitnan'), ...
           min(k_data,[],'all','omitnan'), ...
           max(k_data,[],'all','omitnan'), ...
           std(k_data,0,'all','omitnan'))
   
   figure()
   z = squeeze(mean(z_data,1));
   k = squeeze(mean(k_data,1));
   loglog(abs(lambdas),z,abs(lambdas),k)
   legend('(z)','(k)')
   xlabel('Stiffness: |\lambda|'), ylabel('accuracy')
   title(sprintf('Accuracy vs stiffness (method = %s)',mname))
   saveas(gcf, sprintf('timedep_M-accuracy_vs_stiffness-DIRK_b%i.png',ib))
   
   figure()
   z = squeeze(mean(z_data,2));
   k = squeeze(mean(k_data,2));
   loglog(abs(lambdas),z,abs(lambdas),k)
   legend('(z)','(k)')
   xlabel('Mass units: \gamma'), ylabel('accuracy')
   title(sprintf('Accuracy vs mass units (method = %s)',mname))
   saveas(gcf, sprintf('timedep_M-accuracy_vs_massunits-DIRK_b%i.png',ib))
end

for ib = 1:length(ARKImethods)
   mname1 = ARKEmethods{ib};
   Be = butcher(mname1);  s = numel(Be(1,:))-1;  q = Be(s+1,1);
   mname2 = ARKImethods{ib};
   Bi = butcher(mname2);
   fprintf('\n  ARK integrator: %s/%s (order = %i)\n',mname1,mname2,q)
   z_data = ARK_rmserrs(:,:,ib,1,end);
   k_data = ARK_rmserrs(:,:,ib,2,end);
   fprintf('      (z): mean = %g, min = %g, max = %g, std = %g\n', ...
           mean(z_data,'all','omitnan'), ...
           min(z_data,[],'all','omitnan'), ...
           max(z_data,[],'all','omitnan'), ...
           std(z_data,0,'all','omitnan'))
   fprintf('      (k): mean = %g, min = %g, max = %g, std = %g\n', ...
           mean(k_data,'all','omitnan'), ...
           min(k_data,[],'all','omitnan'), ...
           max(k_data,[],'all','omitnan'), ...
           std(k_data,0,'all','omitnan'))
   
   figure()
   z = squeeze(mean(z_data,1));
   k = squeeze(mean(k_data,1));
   loglog(abs(lambdas),z,abs(lambdas),k)
   legend('(z)','(k)')
   xlabel('Stiffness: |\lambda|'), ylabel('accuracy')
   title(sprintf('Accuracy vs stiffness (method = %s)',mname1))
   saveas(gcf, sprintf('timedep_M-accuracy_vs_stiffness-ARK_b%i.png',ib))
   
   figure()
   z = squeeze(mean(z_data,2));
   k = squeeze(mean(k_data,2));
   loglog(abs(lambdas),z,abs(lambdas),k)
   legend('(z)','(k)')
   xlabel('Mass units: \gamma'), ylabel('accuracy')
   title(sprintf('Accuracy vs mass units (method = %s)',mname1))
   saveas(gcf, sprintf('timedep_M-accuracy_vs_massunits-ARK_b%i.png',ib))
end

nlam = sum(abs(lambdas) <= 100);
for ib = 1:length(ERKmethods)
   mname = ERKmethods{ib};
   B = butcher(mname);  s = numel(B(1,:))-1;  q = B(s+1,1);
   fprintf('\n  ERK integrator: %s (order = %i)\n',mname,q)
   z_data = ERK_rmserrs(:,1:nlam,ib,1,end);
   k_data = ERK_rmserrs(:,1:nlam,ib,2,end);
   fprintf('      (z): mean = %g, min = %g, max = %g, std = %g\n', ...
           mean(z_data,'all','omitnan'), ...
           min(z_data,[],'all','omitnan'), ...
           max(z_data,[],'all','omitnan'), ...
           std(z_data,0,'all','omitnan'))
   fprintf('      (k): mean = %g, min = %g, max = %g, std = %g\n', ...
           mean(k_data,'all','omitnan'), ...
           min(k_data,[],'all','omitnan'), ...
           max(k_data,[],'all','omitnan'), ...
           std(k_data,0,'all','omitnan'))
   
   figure()
   z = squeeze(mean(z_data,1));
   k = squeeze(mean(k_data,1));
   loglog(abs(lambdas(1:nlam)),z,abs(lambdas(1:nlam)),k)
   legend('(z)','(k)')
   xlabel('Stiffness: |\lambda|'), ylabel('accuracy')
   title(sprintf('Accuracy vs stiffness (method = %s)',mname))
   saveas(gcf, sprintf('timedep_M-accuracy_vs_stiffness-ERK_b%i.png',ib))
   
   figure()
   z = squeeze(mean(z_data,2));
   k = squeeze(mean(k_data,2));
   loglog(abs(lambdas),z,abs(lambdas),k)
   legend('(z)','(k)')
   xlabel('Mass units: \gamma'), ylabel('accuracy')
   title(sprintf('Accuracy vs mass units (method = %s)',mname))
   saveas(gcf, sprintf('timedep_M-accuracy_vs_massunits-ERK_b%i.png',ib))
end



% Newton iterations at smallest step size
fprintf('\nNewton iterations at smallest step size:\n');

for ib = 1:length(DIRKmethods)
   mname = DIRKmethods{ib};
   B = butcher(mname);  s = numel(B(1,:))-1;  q = B(s+1,1);
   fprintf('\n  DIRK integrator: %s (order = %i)\n',mname,q)
   z_data = DIRK_lsolves(:,:,ib,1,end);
   k_data = DIRK_lsolves(:,:,ib,2,end);
   fprintf('      (z): mean = %g, min = %g, max = %g, std = %g\n', ...
           mean(z_data,'all','omitnan'), ...
           min(z_data,[],'all','omitnan'), ...
           max(z_data,[],'all','omitnan'), ...
           std(z_data,0,'all','omitnan'))
   fprintf('      (k): mean = %g, min = %g, max = %g, std = %g\n', ...
           mean(k_data,'all','omitnan'), ...
           min(k_data,[],'all','omitnan'), ...
           max(k_data,[],'all','omitnan'), ...
           std(k_data,0,'all','omitnan'))
   
   figure()
   z = squeeze(mean(z_data,1));
   k = squeeze(mean(k_data,1));
   semilogx(abs(lambdas),z,abs(lambdas),k)
   legend('(z)','(k)')
   xlabel('Stiffness: |\lambda|'), ylabel('NIters')
   title(sprintf('Newt Iters vs stiffness (method = %s)',mname))
   saveas(gcf, sprintf('timedep_M-newton_vs_stiffness-DIRK_b%i.png',ib))
   
   figure()
   z = squeeze(mean(z_data,2));
   k = squeeze(mean(k_data,2));
   loglog(abs(lambdas),z,abs(lambdas),k)
   legend('(z)','(k)')
   xlabel('Mass units: \gamma'), ylabel('NIters')
   title(sprintf('Newt Iters vs stiffness (method = %s)',mname))
   saveas(gcf, sprintf('timedep_M-newton_vs_stiffness-DIRK_b%i.png',ib))
end

for ib = 1:length(ARKImethods)
   mname1 = ARKEmethods{ib};
   Be = butcher(mname1);  s = numel(Be(1,:))-1;  q = Be(s+1,1);
   mname2 = ARKImethods{ib};
   Bi = butcher(mname2);
   fprintf('\n  ARK integrator: %s/%s (order = %i)\n',mname1,mname2,q)
   z_data = ARK_lsolves(:,:,ib,1,end);
   k_data = ARK_lsolves(:,:,ib,2,end);
   fprintf('      (z): mean = %g, min = %g, max = %g, std = %g\n', ...
           mean(z_data,'all','omitnan'), ...
           min(z_data,[],'all','omitnan'), ...
           max(z_data,[],'all','omitnan'), ...
           std(z_data,0,'all','omitnan'))
   fprintf('      (k): mean = %g, min = %g, max = %g, std = %g\n', ...
           mean(k_data,'all','omitnan'), ...
           min(k_data,[],'all','omitnan'), ...
           max(k_data,[],'all','omitnan'), ...
           std(k_data,0,'all','omitnan'))
   
   figure()
   z = squeeze(mean(z_data,1));
   k = squeeze(mean(k_data,1));
   semilogx(abs(lambdas),z,abs(lambdas),k)
   legend('(z)','(k)')
   xlabel('Stiffness: |\lambda|'), ylabel('NIters')
   title(sprintf('Newt Iters vs stiffness (method = %s)',mname1))
   saveas(gcf, sprintf('timedep_M-newton_vs_stiffness-ARK_b%i.png',ib))
   
   figure()
   z = squeeze(mean(z_data,2));
   k = squeeze(mean(k_data,2));
   semilogx(abs(lambdas),z,abs(lambdas),k)
   legend('(z)','(k)')
   xlabel('Mass units: \gamma'), ylabel('NIters')
   title(sprintf('Newt Iters vs stiffness (method = %s)',mname1))
   saveas(gcf, sprintf('timedep_M-newton_vs_stiffness-ARK_b%i.png',ib))
end





%--- M=I vs M=M(t) ---%
fprintf('\n\nComparing M=I vs M=M(t) (for fixed gamma = %g):\n', gammas(end));

% Solution order
fprintf('\nOrder of accuracy:\n');

for ib = 1:length(DIRKmethods)
   mname = DIRKmethods{ib};
   B = butcher(mname);  s = numel(B(1,:))-1;  q = B(s+1,1);
   fprintf('\n  DIRK integrator: %s (order = %i)\n',mname,q)
   z_data = squeeze(DIRK_ord(end,:,ib,:));
   k_data = squeeze(DIRK_ord_noM(:,ib,:));
   fprintf('      M(t): mean = %g, min = %g, max = %g, std = %g\n', ...
           mean(z_data,'all','omitnan'), ...
           min(z_data,[],'all','omitnan'), ...
           max(z_data,[],'all','omitnan'), ...
           std(z_data,0,'all','omitnan'))
   fprintf('         I: mean = %g, min = %g, max = %g, std = %g\n', ...
           mean(k_data,'all','omitnan'), ...
           min(k_data,[],'all','omitnan'), ...
           max(k_data,[],'all','omitnan'), ...
           std(k_data,0,'all','omitnan'))
   
   figure()
   z = squeeze(mean(z_data,2));
   k = squeeze(mean(k_data,2));
   semilogx(abs(lambdas),z,abs(lambdas),k)
   legend('M(t)','I')
   xlabel('Stiffness: |\lambda|'), ylabel('order')
   title(sprintf('Measured order vs stiffness (method = %s)',mname))
   saveas(gcf, sprintf('timedep_M-MvsI-order_vs_stiffness-DIRK_b%i.png',ib))
end

for ib = 1:length(ARKImethods)
   mname1 = ARKEmethods{ib};
   Be = butcher(mname1);  s = numel(Be(1,:))-1;  q = Be(s+1,1);
   mname2 = ARKImethods{ib};
   Bi = butcher(mname2);
   fprintf('\n  ARK integrator: %s/%s (order = %i)\n',mname1,mname2,q)
   z_data = squeeze(ARK_ord(end,:,ib,:));
   k_data = squeeze(ARK_ord_noM(:,ib,:));
   fprintf('      M(t): mean = %g, min = %g, max = %g, std = %g\n', ...
           mean(z_data,'all','omitnan'), ...
           min(z_data,[],'all','omitnan'), ...
           max(z_data,[],'all','omitnan'), ...
           std(z_data,0,'all','omitnan'))
   fprintf('         I: mean = %g, min = %g, max = %g, std = %g\n', ...
           mean(k_data,'all','omitnan'), ...
           min(k_data,[],'all','omitnan'), ...
           max(k_data,[],'all','omitnan'), ...
           std(k_data,0,'all','omitnan'))
   
   figure()
   z = squeeze(mean(z_data,2));
   k = squeeze(mean(k_data,2));
   semilogx(abs(lambdas),z,abs(lambdas),k)
   legend('M(t)','I')
   xlabel('Stiffness: |\lambda|'), ylabel('order')
   title(sprintf('Measured order vs stiffness (method = %s)',mname1))
   saveas(gcf, sprintf('timedep_M-MvsI-order_vs_stiffness-ARK_b%i.png',ib))
end

nlam = sum(abs(lambdas) <= 100);
for ib = 1:length(ERKmethods)
   mname = ERKmethods{ib};
   B = butcher(mname);  s = numel(B(1,:))-1;  q = B(s+1,1);
   fprintf('\n  ERK integrator: %s (order = %i)\n',mname,q)
   z_data = squeeze(ERK_ord(end,1:nlam,ib,:));
   k_data = squeeze(ERK_ord_noM(1:nlam,ib,:));
   fprintf('      M(t): mean = %g, min = %g, max = %g, std = %g\n', ...
           mean(z_data,'all','omitnan'), ...
           min(z_data,[],'all','omitnan'), ...
           max(z_data,[],'all','omitnan'), ...
           std(z_data,0,'all','omitnan'))
   fprintf('         I: mean = %g, min = %g, max = %g, std = %g\n', ...
           mean(k_data,'all','omitnan'), ...
           min(k_data,[],'all','omitnan'), ...
           max(k_data,[],'all','omitnan'), ...
           std(k_data,0,'all','omitnan'))
   
   figure()
   z = squeeze(mean(z_data,2));
   k = squeeze(mean(k_data,2));
   semilogx(abs(lambdas(1:nlam)),z,abs(lambdas(1:nlam)),k)
   legend('M(t)','I')
   xlabel('Stiffness: |\lambda|'), ylabel('order')
   title(sprintf('Measured order vs stiffness (method = %s)',mname))
   saveas(gcf, sprintf('timedep_M-MvsI-order_vs_stiffness-ERK_b%i.png',ib))
end



% Solution accuracy
fprintf('\nAccuracy at smallest step size:\n');

for ib = 1:length(DIRKmethods)
   mname = DIRKmethods{ib};
   B = butcher(mname);  s = numel(B(1,:))-1;  q = B(s+1,1);
   fprintf('\n  DIRK integrator: %s (order = %i)\n',mname,q)
   z_data = squeeze(DIRK_rmserrs(end,:,ib,:,end));
   k_data = squeeze(DIRK_rmserrs_noM(:,ib,:,end));
   fprintf('      M(t): mean = %g, min = %g, max = %g, std = %g\n', ...
           mean(z_data,'all','omitnan'), ...
           min(z_data,[],'all','omitnan'), ...
           max(z_data,[],'all','omitnan'), ...
           std(z_data,0,'all','omitnan'))
   fprintf('         I: mean = %g, min = %g, max = %g, std = %g\n', ...
           mean(k_data,'all','omitnan'), ...
           min(k_data,[],'all','omitnan'), ...
           max(k_data,[],'all','omitnan'), ...
           std(k_data,0,'all','omitnan'))
   
   figure()
   z = squeeze(mean(z_data,2));
   k = squeeze(mean(k_data,2));
   loglog(abs(lambdas),z,abs(lambdas),k)
   legend('M(t)','I')
   xlabel('Stiffness: |\lambda|'), ylabel('accuracy')
   title(sprintf('Accuracy vs stiffness (method = %s)',mname))
   saveas(gcf, sprintf('timedep_M-MvsI-accuracy_vs_stiffness-DIRK_b%i.png',ib))
end

for ib = 1:length(ARKImethods)
   mname1 = ARKEmethods{ib};
   Be = butcher(mname1);  s = numel(Be(1,:))-1;  q = Be(s+1,1);
   mname2 = ARKImethods{ib};
   Bi = butcher(mname2);
   fprintf('\n  ARK integrator: %s/%s (order = %i)\n',mname1,mname2,q)
   z_data = squeeze(ARK_rmserrs(end,:,ib,:,end));
   k_data = squeeze(ARK_rmserrs_noM(:,ib,:,end));
   fprintf('      M(t): mean = %g, min = %g, max = %g, std = %g\n', ...
           mean(z_data,'all','omitnan'), ...
           min(z_data,[],'all','omitnan'), ...
           max(z_data,[],'all','omitnan'), ...
           std(z_data,0,'all','omitnan'))
   fprintf('         I: mean = %g, min = %g, max = %g, std = %g\n', ...
           mean(k_data,'all','omitnan'), ...
           min(k_data,[],'all','omitnan'), ...
           max(k_data,[],'all','omitnan'), ...
           std(k_data,0,'all','omitnan'))
   
   figure()
   z = squeeze(mean(z_data,2));
   k = squeeze(mean(k_data,2));
   loglog(abs(lambdas),z,abs(lambdas),k)
   legend('M(t)','I')
   xlabel('Stiffness: |\lambda|'), ylabel('accuracy')
   title(sprintf('Accuracy vs stiffness (method = %s)',mname1))
   saveas(gcf, sprintf('timedep_M-MvsI-accuracy_vs_stiffness-ARK_b%i.png',ib))
end

nlam = sum(abs(lambdas) <= 100);
for ib = 1:length(ERKmethods)
   mname = ERKmethods{ib};
   B = butcher(mname);  s = numel(B(1,:))-1;  q = B(s+1,1);
   fprintf('\n  ERK integrator: %s (order = %i)\n',mname,q)
   z_data = squeeze(ERK_rmserrs(end,1:nlam,ib,:,end));
   k_data = squeeze(ERK_rmserrs_noM(1:nlam,ib,:,end));
   fprintf('      M(t): mean = %g, min = %g, max = %g, std = %g\n', ...
           mean(z_data,'all','omitnan'), ...
           min(z_data,[],'all','omitnan'), ...
           max(z_data,[],'all','omitnan'), ...
           std(z_data,0,'all','omitnan'))
   fprintf('         I: mean = %g, min = %g, max = %g, std = %g\n', ...
           mean(k_data,'all','omitnan'), ...
           min(k_data,[],'all','omitnan'), ...
           max(k_data,[],'all','omitnan'), ...
           std(k_data,0,'all','omitnan'))
   
   figure()
   z = squeeze(mean(z_data,2));
   k = squeeze(mean(k_data,2));
   loglog(abs(lambdas(1:nlam)),z,abs(lambdas(1:nlam)),k)
   legend('M(t)','I')
   xlabel('Stiffness: |\lambda|'), ylabel('accuracy')
   title(sprintf('Accuracy vs stiffness (method = %s)',mname))
   saveas(gcf, sprintf('timedep_M-MvsI-accuracy_vs_stiffness-ERK_b%i.png',ib))
end



% Newton iterations at smallest step size
fprintf('\nNewton iterations at smallest step size:\n');

for ib = 1:length(DIRKmethods)
   mname = DIRKmethods{ib};
   B = butcher(mname);  s = numel(B(1,:))-1;  q = B(s+1,1);
   fprintf('\n  DIRK integrator: %s (order = %i)\n',mname,q)
   z_data = squeeze(DIRK_lsolves(end,:,ib,:,end));
   k_data = squeeze(DIRK_lsolves_noM(:,ib,:,end));
   fprintf('      M(t): mean = %g, min = %g, max = %g, std = %g\n', ...
           mean(z_data,'all','omitnan'), ...
           min(z_data,[],'all','omitnan'), ...
           max(z_data,[],'all','omitnan'), ...
           std(z_data,0,'all','omitnan'))
   fprintf('         I: mean = %g, min = %g, max = %g, std = %g\n', ...
           mean(k_data,'all','omitnan'), ...
           min(k_data,[],'all','omitnan'), ...
           max(k_data,[],'all','omitnan'), ...
           std(k_data,0,'all','omitnan'))
   
   figure()
   z = squeeze(mean(z_data,2));
   k = squeeze(mean(k_data,2));
   semilogx(abs(lambdas),z,abs(lambdas),k)
   legend('M(t)','I')
   xlabel('Stiffness: |\lambda|'), ylabel('NIters')
   title(sprintf('Newt Iters vs stiffness (method = %s)',mname))
   saveas(gcf, sprintf('timedep_M-MvsI-newton_vs_stiffness-DIRK_b%i.png',ib))
end

for ib = 1:length(ARKImethods)
   mname1 = ARKEmethods{ib};
   Be = butcher(mname1);  s = numel(Be(1,:))-1;  q = Be(s+1,1);
   mname2 = ARKImethods{ib};
   Bi = butcher(mname2);
   fprintf('\n  ARK integrator: %s/%s (order = %i)\n',mname1,mname2,q)
   z_data = squeeze(ARK_lsolves(end,:,ib,:,end));
   k_data = squeeze(ARK_lsolves_noM(:,ib,:,end));
   fprintf('      M(t): mean = %g, min = %g, max = %g, std = %g\n', ...
           mean(z_data,'all','omitnan'), ...
           min(z_data,[],'all','omitnan'), ...
           max(z_data,[],'all','omitnan'), ...
           std(z_data,0,'all','omitnan'))
   fprintf('         I: mean = %g, min = %g, max = %g, std = %g\n', ...
           mean(k_data,'all','omitnan'), ...
           min(k_data,[],'all','omitnan'), ...
           max(k_data,[],'all','omitnan'), ...
           std(k_data,0,'all','omitnan'))
   
   figure()
   z = squeeze(mean(z_data,2));
   k = squeeze(mean(k_data,2));
   semilogx(abs(lambdas),z,abs(lambdas),k)
   legend('M(t)','I')
   xlabel('Stiffness: |\lambda|'), ylabel('NIters')
   title(sprintf('Newt Iters vs stiffness (method = %s)',mname1))
   saveas(gcf, sprintf('timedep_M-MvsI-newton_vs_stiffness-ARK_b%i.png',ib))
end

% close diary
diary off

% end of script
