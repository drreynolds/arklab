% Driver to analyze output from the script 'driver_timedep_bdry.m'.
% The primary purpose of that script, and these analyses, is to
% determine whether there is any appreciable performance difference
% between approaches that enforce time-dependent boundary
% conditions through the ODE RHS function ('recommended'), or by 
% forcing the values directly into each stage ('forced').  We thus 
% examine a variety of metrics to assess the performance of both
% methods.
%
% Daniel R. Reynolds
% Department of Mathematics
% Southern Methodist University
% November 2019
% All Rights Reserved
clear

% start output diary
!\rm timedep_bdry_analysis.txt
diary timedep_bdry_analysis.txt

% load results:
%   Test setup:
%     lambdas = stiffness parameters
%     DIRKmethods = DIRK Butcher tables (strings)
%     ARKEmethods = ARK explicit Butcher tables (strings)
%     ARKImethods = ARK ipmplicit Butcher tables (strings)
%     ERKmethods = ERK Butcher tables (strings)
%   'Recommended' Results:
%     DIRK_recommended_rmserrs(lambdas,DIRKmethods,hvals)
%     DIRK_recommended_ord(lambdas,DIRKmethods)
%     DIRK_recommended_ordred(lambdas,DIRKmethods)
%     ARK_recommended_rmserrs(lambdas,ARKmethods,hvals)
%     ARK_recommended_ord(lambdas,ARKmethods)
%     ARK_recommended_ordred(lambdas,ARKmethods)
%     ERK_recommended_rmserrs(lambdas,ERKmethods,hvals)
%     ERK_recommended_ord(lambdas,ERKmethods)
%     ERK_recommended_ordred(lambdas,ERKmethods)
%   'Forced' Results:
%     DIRK_forced_rmserrs(lambdas,DIRKmethods,hvals)
%     DIRK_forced_ord(lambdas,DIRKmethods)
%     DIRK_forced_ordred(lambdas,DIRKmethods)
%     ARK_forced_rmserrs(lambdas,ARKmethods,hvals)
%     ARK_forced_ord(lambdas,ARKmethods)
%     ARK_forced_ordred(lambdas,ARKmethods)
%     ERK_forced_rmserrs(lambdas,ERKmethods,hvals)
%     ERK_forced_ord(lambdas,ERKmethods)
%     ERK_forced_ordred(lambdas,ERKmethods)
load timedep_bdry_data.mat

% set shortcut variables for testing array lengths
nl = length(lambdas);
nh = length(hvals);
nDIRK = length(DIRKmethods);
nARK = length(ARKEmethods);
nERK = length(ERKmethods);

fprintf('\n\nComparing recommended vs forced approaches for time-dependent boundary conditions:\n');


% Solution order
fprintf('\nOrder of accuracy:\n');

for ib = 1:length(DIRKmethods)
   mname = DIRKmethods{ib};
   B = butcher(mname);  s = numel(B(1,:))-1;  q = B(s+1,1);
   fprintf('\n  DIRK integrator: %s (order = %i)\n',mname,q)
   r_data = DIRK_recommended_ord(:,ib);
   f_data = DIRK_forced_ord(:,ib);
   fprintf('      recommended: mean = %g, min = %g, max = %g, std = %g\n', ...
           mean(r_data,'omitnan'), ...
           min(r_data,[],'omitnan'), ...
           max(r_data,[],'omitnan'), ...
           std(r_data,'omitnan'))
   fprintf('      forced: mean = %g, min = %g, max = %g, std = %g\n', ...
           mean(f_data,'omitnan'), ...
           min(f_data,[],'omitnan'), ...
           max(f_data,[],'omitnan'), ...
           std(f_data,'omitnan'))
   
   figure()
   semilogx(abs(lambdas),r_data,abs(lambdas),f_data)
   legend('recommended','forced')
   xlabel('Stiffness: |\lambda|'), ylabel('order')
   title(sprintf('Measured order vs stiffness (method = %s)',mname))
   saveas(gcf, sprintf('timedep_bdry-order_vs_stiffness-DIRK_b%i.png',ib))
end

for ib = 1:length(ARKImethods)
   mname1 = ARKEmethods{ib};
   Be = butcher(mname1);  s = numel(Be(1,:))-1;  q = Be(s+1,1);
   mname2 = ARKImethods{ib};
   Bi = butcher(mname2);
   fprintf('\n    ARK integrator: %s/%s (order = %i)\n',mname1,mname2,q)
   r_data = ARK_recommended_ord(:,ib);
   f_data = ARK_forced_ord(:,ib);
   fprintf('      recommended: mean = %g, min = %g, max = %g, std = %g\n', ...
           mean(r_data,'omitnan'), ...
           min(r_data,[],'omitnan'), ...
           max(r_data,[],'omitnan'), ...
           std(r_data,'omitnan'))
   fprintf('      forced: mean = %g, min = %g, max = %g, std = %g\n', ...
           mean(f_data,'omitnan'), ...
           min(f_data,[],'omitnan'), ...
           max(f_data,[],'omitnan'), ...
           std(f_data,'omitnan'))
   
   figure()
   semilogx(abs(lambdas),r_data,abs(lambdas),f_data)
   legend('recommended','forced')
   xlabel('Stiffness: |\lambda|'), ylabel('order')
   title(sprintf('Measured order vs stiffness (method = %s)',mname1))
   saveas(gcf, sprintf('timedep_bdry-order_vs_stiffness-ARK_b%i.png',ib))
end

nlam = sum(abs(lambdas) <= 100);
for ib = 1:length(ERKmethods)
   mname = ERKmethods{ib};
   B = butcher(mname);  s = numel(B(1,:))-1;  q = B(s+1,1);
   fprintf('\n  ERK integrator: %s (order = %i)\n',mname,q)
   r_data = ERK_recommended_ord(1:nlam,ib);
   f_data = ERK_forced_ord(1:nlam,ib);
   fprintf('      recommended: mean = %g, min = %g, max = %g, std = %g\n', ...
           mean(r_data,'omitnan'), ...
           min(r_data,[],'omitnan'), ...
           max(r_data,[],'omitnan'), ...
           std(r_data,'omitnan'))
   fprintf('      forced: mean = %g, min = %g, max = %g, std = %g\n', ...
           mean(f_data,'omitnan'), ...
           min(f_data,[],'omitnan'), ...
           max(f_data,[],'omitnan'), ...
           std(f_data,'omitnan'))
   
   figure()
   semilogx(abs(lambdas(1:nlam)),r_data,abs(lambdas(1:nlam)),f_data)
   legend('recommended','forced')
   xlabel('Stiffness: |\lambda|'), ylabel('order')
   title(sprintf('Measured order vs stiffness (method = %s)',mname))
   saveas(gcf, sprintf('timedep_bdry-order_vs_stiffness-ERK_b%i.png',ib))
end




% Solution accuracy
fprintf('\nAccuracy at smallest step size:\n');

for ib = 1:length(DIRKmethods)
   mname = DIRKmethods{ib};
   B = butcher(mname);  s = numel(B(1,:))-1;  q = B(s+1,1);
   fprintf('\n  DIRK integrator: %s (order = %i)\n',mname,q)
   r_data = DIRK_recommended_rmserrs(:,ib,end);
   f_data = DIRK_forced_rmserrs(:,ib,end);
   fprintf('      recommended: mean = %g, min = %g, max = %g, std = %g\n', ...
           mean(r_data,'omitnan'), ...
           min(r_data,[],'omitnan'), ...
           max(r_data,[],'omitnan'), ...
           std(r_data,'omitnan'))
   fprintf('      forced: mean = %g, min = %g, max = %g, std = %g\n', ...
           mean(f_data,'omitnan'), ...
           min(f_data,[],'omitnan'), ...
           max(f_data,[],'omitnan'), ...
           std(f_data,'omitnan'))
   
   figure()
   loglog(abs(lambdas),r_data,abs(lambdas),f_data)
   legend('recommended','forced')
   xlabel('Stiffness: |\lambda|'), ylabel('accuracy')
   title(sprintf('Accuracy vs stiffness (method = %s)',mname))
   saveas(gcf, sprintf('timedep_bdry-accuracy_vs_stiffness-DIRK_b%i.png',ib))
end

for ib = 1:length(ARKImethods)
   mname1 = ARKEmethods{ib};
   Be = butcher(mname1);  s = numel(Be(1,:))-1;  q = Be(s+1,1);
   mname2 = ARKImethods{ib};
   Bi = butcher(mname2);
   fprintf('\n    ARK integrator: %s/%s (order = %i)\n',mname1,mname2,q)
   r_data = ARK_recommended_rmserrs(:,ib,end);
   f_data = ARK_forced_rmserrs(:,ib,end);
   fprintf('      recommended: mean = %g, min = %g, max = %g, std = %g\n', ...
           mean(r_data,'omitnan'), ...
           min(r_data,[],'omitnan'), ...
           max(r_data,[],'omitnan'), ...
           std(r_data,'omitnan'))
   fprintf('      forced: mean = %g, min = %g, max = %g, std = %g\n', ...
           mean(f_data,'omitnan'), ...
           min(f_data,[],'omitnan'), ...
           max(f_data,[],'omitnan'), ...
           std(f_data,'omitnan'))
   
   figure()
   loglog(abs(lambdas),r_data,abs(lambdas),f_data)
   legend('recommended','forced')
   xlabel('Stiffness: |\lambda|'), ylabel('accuracy')
   title(sprintf('Accuracy vs stiffness (method = %s)',mname1))
   saveas(gcf, sprintf('timedep_bdry-accuracy_vs_stiffness-ARK_b%i.png',ib))
end

nlam = sum(abs(lambdas) <= 100);
for ib = 1:length(ERKmethods)
   mname = ERKmethods{ib};
   B = butcher(mname);  s = numel(B(1,:))-1;  q = B(s+1,1);
   fprintf('\n  ERK integrator: %s (order = %i)\n',mname,q)
   r_data = ERK_recommended_rmserrs(1:nlam,ib,end);
   f_data = ERK_forced_rmserrs(1:nlam,ib,end);
   fprintf('      recommended: mean = %g, min = %g, max = %g, std = %g\n', ...
           mean(r_data,'omitnan'), ...
           min(r_data,[],'omitnan'), ...
           max(r_data,[],'omitnan'), ...
           std(r_data,'omitnan'))
   fprintf('      forced: mean = %g, min = %g, max = %g, std = %g\n', ...
           mean(f_data,'omitnan'), ...
           min(f_data,[],'omitnan'), ...
           max(f_data,[],'omitnan'), ...
           std(f_data,'omitnan'))
   
   figure()
   loglog(abs(lambdas(1:nlam)),r_data,abs(lambdas(1:nlam)),f_data)
   legend('recommended','forced')
   xlabel('Stiffness: |\lambda|'), ylabel('accuracy')
   title(sprintf('Accuracy vs stiffness (method = %s)',mname))
   saveas(gcf, sprintf('timedep_bdry-accuracy_vs_stiffness-ERK_b%i.png',ib))
end

% close diary
diary off

% end of script
