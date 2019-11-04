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

% load results:
%   Test setup:
%     lambdas = stiffness parameters
%     gammas = mass matrix units
%     hvals = time step sizes
%     algs = solve for z (0) vs solve for k (1)
%     DIRKmethods = DIRK Butcher tables (strings)
%     ARKEmethods = ARK explicit Butcher tables (strings)
%     ARKImethods = ARK ipmplicit Butcher tables (strings)
%     ERKmethods = ERK Butcher tables (strings)
%   Results for M=M(t):
%     DIRK_rmserrs(gammas,lambdas,DIRKmethods,algs,hvals)
%     DIRK_lsolves(gammas,lambdas,DIRKmethods,algs,hvals)
%     DIRK_ord(gammas,lambdas,DIRKmethods,algs,hvals)
%     DIRK_ordred(gammas,lambdas,DIRKmethods,algs,hvals)
%     ARK_rmserrs(gammas,lambdas,DIRKmethods,algs,hvals)
%     ARK_lsolves(gammas,lambdas,DIRKmethods,algs,hvals)
%     ARK_ord(gammas,lambdas,DIRKmethods,algs,hvals)
%     ARK_ordred(gammas,lambdas,DIRKmethods,algs,hvals)
%     ERK_rmserrs(gammas,lambdas,DIRKmethods,algs,hvals)
%     ERK_ord(gammas,lambdas,DIRKmethods,algs,hvals)
%     ERK_ordred(gammas,lambdas,DIRKmethods,algs,hvals)
%   Results for M=I:
%     DIRK_rmserrs_noM(lambdas,DIRKmethods,algs,hvals)
%     DIRK_lsolves_noM(lambdas,DIRKmethods,algs,hvals)
%     DIRK_ord_noM(lambdas,DIRKmethods,algs,hvals)
%     DIRK_ordred_noM(lambdas,DIRKmethods,algs,hvals)
%     ARK_rmserrs_noM(lambdas,DIRKmethods,algs,hvals)
%     ARK_lsolves_noM(lambdas,DIRKmethods,algs,hvals)
%     ARK_ord_noM(lambdas,DIRKmethods,algs,hvals)
%     ARK_ordred_noM(lambdas,DIRKmethods,algs,hvals)
%     ERK_rmserrs_noM(lambdas,DIRKmethods,algs,hvals)
%     ERK_ord_noM(lambdas,DIRKmethods,algs,hvals)
%     ERK_ordred_noM(lambdas,DIRKmethods,algs,hvals)
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

% Solution accuracy
fprintf('\n  Solution accuracy ratios (z)./(k):\n');

%    DIRK solver results
d_data = DIRK_rmserrs(:,:,:,1,:)./DIRK_rmserrs(:,:,:,2,:);
fprintf('    DIRK: mean = %g, min = %g, max = %g, std = %g\n', ...
        mean(d_data,'all','omitnan'), ...
        min(d_data,[],'all','omitnan'), ...
        max(d_data,[],'all','omitnan'), ...
        std(d_data,0,'all','omitnan'))

%    ARK solver results
a_data = ARK_rmserrs(:,:,:,1,:)./ARK_rmserrs(:,:,:,2,:);
fprintf('    ARK:  mean = %g, min = %g, max = %g, std = %g\n', ...
        mean(a_data,'all','omitnan'), ...
        min(a_data,[],'all','omitnan'), ...
        max(a_data,[],'all','omitnan'), ...
        std(a_data,0,'all'))

%    ERK solver results
e_data = ERK_rmserrs(:,:,:,1,:)./ERK_rmserrs(:,:,:,2,:);
fprintf('    ERK:  mean = %g, min = %g, max = %g, std = %g\n', ...
        mean(e_data,'all','omitnan'), ...
        min(e_data,[],'all','omitnan'), ...
        max(e_data,[],'all','omitnan'), ...
        std(e_data,0,'all','omitnan'))

%    accuracy vs stiffness (average over all other parameters)
figure()
d = squeeze(mean(mean(mean(d_data,5),3),1));
a = squeeze(mean(mean(mean(a_data,5),3),1));
e = squeeze(mean(mean(mean(e_data,5),3),1));
semilogx(abs(lambdas),d,abs(lambdas),a,abs(lambdas(1:length(e))),e)
legend('DIRK','ARK','ERK')
xlabel('Stiffness: |\lambda|'), ylabel('(z)./(k)')
title('Accuracy ratio vs stiffness')
saveas(gcf, 'z_vs_k-error_vs_stiffness.png')

%    accuracy vs mass matrix units (average over all other parameters)
figure()
d = squeeze(mean(mean(mean(d_data,5),3),2));
a = squeeze(mean(mean(mean(a_data,5),3),2));
e = squeeze(mean(mean(mean(e_data,5),3),2));
semilogx(gammas,d,gammas,a,gammas,e)
legend('DIRK','ARK','ERK')
xlabel('Mass units: \gamma'), ylabel('(z)./(k)')
title('Accuracy ratio vs mass matrix units')
saveas(gcf, 'z_vs_k-error_vs_massunits.png')


% Newton iterations
fprintf('\n  Newton iteration ratios (z)./(k):\n');

%    DIRK solver results
d_data = DIRK_lsolves(:,:,:,1,:)./DIRK_lsolves(:,:,:,2,:);
fprintf('    DIRK: mean = %g, min = %g, max = %g, std = %g\n', ...
        mean(d_data,'all','omitnan'), ...
        min(d_data,[],'all','omitnan'), ...
        max(d_data,[],'all','omitnan'), ...
        std(d_data,0,'all','omitnan'))

%    ARK solver results
a_data = ARK_lsolves(:,:,:,1,:)./ARK_lsolves(:,:,:,2,:);
fprintf('    ARK:  mean = %g, min = %g, max = %g, std = %g\n', ...
        mean(a_data,'all','omitnan'), ...
        min(a_data,[],'all','omitnan'), ...
        max(a_data,[],'all','omitnan'), ...
        std(a_data,0,'all','omitnan'))

%    Newton iterations vs stiffness (average over all other parameters)
figure()
d = squeeze(mean(mean(mean(d_data,5),3),1));
a = squeeze(mean(mean(mean(a_data,5),3),1));
semilogx(abs(lambdas),d,abs(lambdas),a)
legend('DIRK','ARK')
xlabel('Stiffness: |\lambda|'), ylabel('(z)./(k)')
title('Newton iteration ratio vs stiffness')
saveas(gcf, 'z_vs_k-newton_vs_stiffness.png')

%    Newton iterations vs mass matrix units (average over all other parameters)
figure()
d = squeeze(mean(mean(mean(d_data,5),3),2));
a = squeeze(mean(mean(mean(a_data,5),3),2));
semilogx(gammas,d,gammas,a)
legend('DIRK','ARK')
xlabel('Mass units: \gamma'), ylabel('(z)./(k)')
title('Newton iteration ratio vs mass matrix units')
saveas(gcf, 'z_vs_k-newton_vs_massunits.png')


% order reduction
fprintf('\n  Order reduction difference (z) - (k):\n');

%    DIRK solver results
d_data = DIRK_ordred(:,:,:,1,:) - DIRK_ordred(:,:,:,2,:);
fprintf('    DIRK: mean = %g, min = %g, max = %g, std = %g\n', ...
        mean(d_data,'all','omitnan'), ...
        min(d_data,[],'all','omitnan'), ...
        max(d_data,[],'all','omitnan'), ...
        std(d_data,0,'all','omitnan'))

%    ARK solver results
a_data = ARK_ordred(:,:,:,1,:) - ARK_ordred(:,:,:,2,:);
fprintf('    ARK:  mean = %g, min = %g, max = %g, std = %g\n', ...
        mean(a_data,'all','omitnan'), ...
        min(a_data,[],'all','omitnan'), ...
        max(a_data,[],'all','omitnan'), ...
        std(a_data,0,'all','omitnan'))

%    ERK solver results
e_data = ERK_ordred(:,:,:,1,:) - ERK_ordred(:,:,:,2,:);
fprintf('    ERK:  mean = %g, min = %g, max = %g, std = %g\n', ...
        mean(e_data,'all','omitnan'), ...
        min(e_data,[],'all','omitnan'), ...
        max(e_data,[],'all','omitnan'), ...
        std(e_data,0,'all','omitnan'))

%    order reduction vs stiffness (average over all other parameters)
figure()
d = squeeze(mean(mean(mean(d_data,5),3),1));
a = squeeze(mean(mean(mean(a_data,5),3),1));
e = squeeze(mean(mean(mean(e_data,5),3),1));
semilogx(abs(lambdas),d,abs(lambdas),a,abs(lambdas(1:length(e))),e)
legend('DIRK','ARK','ERK')
xlabel('Stiffness: |\lambda|'), ylabel('(z) - (k)')
title('Order reduction difference vs stiffness')
saveas(gcf, 'z_vs_k-ordred_vs_stiffness.png')

%    order reduction vs mass matrix units (average over all other parameters)
figure()
d = squeeze(mean(mean(mean(d_data,5),3),2));
a = squeeze(mean(mean(mean(a_data,5),3),2));
e = squeeze(mean(mean(mean(e_data,5),3),2));
semilogx(gammas,d,gammas,a,gammas,e)
legend('DIRK','ARK','ERK')
xlabel('Mass units: \gamma'), ylabel('(z) - (k)')
title('Order reduction difference vs mass matrix units')
saveas(gcf, 'z_vs_k-ordred_vs_massunits.png')



%--- M=I vs M=M(t) ---%
fprintf('\n\nComparing M=I vs M=M(t) (for fixed gamma = %g):\n', gammas(end));

% Solution accuracy
fprintf('\n  Solution accuracy ratios (M(t))./(I):\n');

%    DIRK solver results
d_data = squeeze(DIRK_rmserrs(end,:,:,:,:))./DIRK_rmserrs_noM(:,:,:,:);
fprintf('    DIRK: mean = %g, min = %g, max = %g, std = %g\n', ...
        mean(d_data,'all','omitnan'), ...
        min(d_data,[],'all','omitnan'), ...
        max(d_data,[],'all','omitnan'), ...
        std(d_data,0,'all','omitnan'))

%    ARK solver results
a_data = squeeze(ARK_rmserrs(end,:,:,:,:))./ARK_rmserrs_noM(:,:,:,:);
fprintf('    ARK:  mean = %g, min = %g, max = %g, std = %g\n', ...
        mean(a_data,'all','omitnan'), ...
        min(a_data,[],'all','omitnan'), ...
        max(a_data,[],'all','omitnan'), ...
        std(a_data,0,'all','omitnan'))

%    ERK solver results
e_data = squeeze(ERK_rmserrs(end,:,:,:,:))./ERK_rmserrs_noM(:,:,:,:);
fprintf('    ERK:  mean = %g, min = %g, max = %g, std = %g\n', ...
        mean(e_data,'all','omitnan'), ...
        min(e_data,[],'all','omitnan'), ...
        max(e_data,[],'all','omitnan'), ...
        std(e_data,0,'all','omitnan'))

%    accuracy vs stiffness (average over all other parameters)
figure()
d = squeeze(mean(mean(mean(d_data,4),3),2));
a = squeeze(mean(mean(mean(a_data,4),3),2));
e = squeeze(mean(mean(mean(e_data,4),3),2));
semilogx(abs(lambdas),d,abs(lambdas),a,abs(lambdas(1:length(e))),e)
legend('DIRK','ARK','ERK')
xlabel('Stiffness: |\lambda|'), ylabel('(M(t))./(I)')
title('Accuracy ratio vs stiffness')
saveas(gcf, 'Mt_vs_I-error_vs_stiffness.png')


% Newton iterations
fprintf('\n  Newton iteration ratios (M(t))./(I):\n');

%    DIRK solver results
d_data = squeeze(DIRK_lsolves(end,:,:,:,:))./DIRK_lsolves_noM(:,:,:,:);
fprintf('    DIRK: mean = %g, min = %g, max = %g, std = %g\n', ...
        mean(d_data,'all','omitnan'), ...
        min(d_data,[],'all','omitnan'), ...
        max(d_data,[],'all','omitnan'), ...
        std(d_data,0,'all','omitnan'))

%    ARK solver results
a_data = squeeze(ARK_lsolves(end,:,:,:,:))./ARK_lsolves_noM(:,:,:,:);
fprintf('    ARK:  mean = %g, min = %g, max = %g, std = %g\n', ...
        mean(a_data,'all','omitnan'), ...
        min(a_data,[],'all','omitnan'), ...
        max(a_data,[],'all','omitnan'), ...
        std(a_data,0,'all','omitnan'))

%    Newton iterations vs stiffness (average over all other parameters)
figure()
d = squeeze(mean(mean(mean(d_data,4),3),2));
a = squeeze(mean(mean(mean(a_data,4),3),2));
semilogx(abs(lambdas),d,abs(lambdas),a)
legend('DIRK','ARK')
xlabel('Stiffness: |\lambda|'), ylabel('(M(t))./(I)')
title('Newton iteration ratio vs stiffness')
saveas(gcf, 'Mt_vs_I-newton_vs_stiffness.png')


% order reduction
fprintf('\n  Order reduction difference (M(t)) - (I):\n');

%    DIRK solver results
d_data = squeeze(DIRK_ordred(end,:,:,:,:)) - DIRK_ordred_noM(:,:,:,:);
fprintf('    DIRK: mean = %g, min = %g, max = %g, std = %g\n', ...
        mean(d_data,'all','omitnan'), ...
        min(d_data,[],'all','omitnan'), ...
        max(d_data,[],'all','omitnan'), ...
        std(d_data,0,'all','omitnan'))

%    ARK solver results
a_data = squeeze(ARK_ordred(end,:,:,:,:)) - ARK_ordred_noM(:,:,:,:);
fprintf('    ARK:  mean = %g, min = %g, max = %g, std = %g\n', ...
        mean(a_data,'all','omitnan'), ...
        min(a_data,[],'all','omitnan'), ...
        max(a_data,[],'all','omitnan'), ...
        std(a_data,0,'all','omitnan'))

%    ERK solver results
e_data = squeeze(ERK_ordred(end,:,:,:,:)) - ERK_ordred_noM(:,:,:,:);
fprintf('    ERK:  mean = %g, min = %g, max = %g, std = %g\n', ...
        mean(e_data,'all','omitnan'), ...
        min(e_data,[],'all','omitnan'), ...
        max(e_data,[],'all','omitnan'), ...
        std(e_data,0,'all','omitnan'))

%    order reduction vs stiffness (average over all other parameters)
figure()
d = squeeze(mean(mean(mean(d_data,4),3),2));
a = squeeze(mean(mean(mean(a_data,4),3),2));
e = squeeze(mean(mean(mean(e_data,4),3),2));
semilogx(abs(lambdas),d,abs(lambdas),a,abs(lambdas(1:length(e))),e)
legend('DIRK','ARK','ERK')
xlabel('Stiffness: |\lambda|'), ylabel('(M(t)) - (I)')
title('Order reduction difference vs stiffness')
saveas(gcf, 'Mt_vs_I-ordred_vs_stiffness.png')

% end of script
