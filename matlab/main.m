addpath('difference_operators');
addpath('helper_functions');

clear

% Solves the free surface problem using SBP-SAT and the FS1 implementation
% stresses are stored on the half-point grid 
% (in FS1 and SBP_STRONG no stresses on the boundary)

% Select scheme, possible choices are:
% FS1        (Tradititional free-surface implementation) 
% SBP_STRONG (strong enforcement of the boundary condition
% SBP_WEAK   (weak enforcement of the boundary condition)
scheme    = 'FS1';

% Select time integrator, possible choices are:
% LF  (Leap-frog, buggy implementation?) 
% RK4 (Low-storage Runge-Kutta 4)
time_integrator = 'RK4';

% Results are written to disk using filename:
filename  = @(scheme,time_integrator,T) sprintf('data/%s6_%s_%d.txt',...
                                                scheme,time_integrator,T);

% Other options
T         = 100;
grids     = [32 64 128 256 512];
alpha     = 1; % penalty parameter when using SBP-WEAK
order     = 6;
a         = 100;
show_plot = false;
CFL       = 0.4;

i = 1;
for n=grids
  fprintf('Running using %d grid points \n',n);
  time = tic;
  [err h_dim] = solver(a,n,T,CFL,show_plot,scheme,time_integrator,order,alpha);
  time = toc(time);
  fprintf('rel. l2-error:  v: %.3f s: %.3f \n',err.v,err.s);
  fprintf('Time to solution: %f \n',time);
  v(i)         = err.v;
  s(i)         = err.s;
  h(i)         = h_dim*a;
  wallclock(i) = time;
  gridnum(i)   = i;
  i = i +1;
end

% Compute convergence rates
rate_v = [0 log2(v(1:end-1)./v(2:end))./log2(h(1:end-1)./h(2:end))];
rate_s = [0 log2(s(1:end-1)./s(2:end))./log2(h(1:end-1)./h(2:end))];

f = fopen(filename(scheme,time_integrator,T),'w');
fprintf(f,'grid h v s rate_v rate_s wallclock \n');
for i=1:length(grids)
  fprintf(f,'%d %g %g %g %g %g %g \n',gridnum(i),h(i),v(i),s(i),rate_v(i),rate_s(i),wallclock(i));
end
fclose(f);

loglog(wallclock,v,'o-');
xlabel('Time to solution (s)');
ylabel('Relative error');
