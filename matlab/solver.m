function [err h] = solver(a,n,T,CFL,show_plot,scheme,time_integrator,order,alpha)

% Solves the free surface problem using SBP-SAT and the FS1 implementation
% stresses are stored on the half-point grid, 
% For FS1:
% stress: n grid points, velocity: n + 1 grid points
% For SBP-SAT:
% stress: n+2 grid points, velocity: n + 1 grid points
% For SBP strong:
% stress: n grid points, velocity: n + 1 grid points

assert(mod(T,2)==0,'Period must be a multiple of two');


% Possible spatial schemes
FS1 = 0; SBP_WEAK = 0; SBP_STRONG = 0;
% Possible temporal schemes
LF = 0; RK4 = 0;
switch scheme
  case 'FS1'
    FS1    = 1;
  case 'SBP_WEAK'
    SBP_WEAK = 1;
  case 'SBP_STRONG'
    SBP_STRONG = 1;
end

switch time_integrator
  case 'LF'
    LF  = 1;
  case  'RK4'
    RK4 = 1;
end


% Wavespeed
rho       = 1;
mu        = 1;
c         = sqrt(mu/rho); 


L     = 1;
test  = true;
h     = L/n;
dt    = CFL*h/c;
Z     = alpha*rho*c;
rhoi  = 1/rho;


% Construct difference operators
if FS1
  xp = linspace(0,L,n+1)';
  xm = h/2 + xp(1:end-1);
  [Dp Dm] = fs1(n,h);
elseif SBP_WEAK
  [xp,xm,Pp,Pm,Qp,Qm] = sbp_weak(order,n,h,false);
elseif SBP_STRONG
  [xp,xm,Pp,Pm,Qp,Qm] = sbp_strong(order,n,h);
end

xp = xp(:); xm = xm(:);

if SBP_WEAK || SBP_STRONG
  Ppi = inv(Pp);
  Pmi = inv(Pm);
  Dp  = inv(Pp)*Qp;
  Dm  = inv(Pm)*Qm;
end

% Construct semi-discrete difference approximation
np   = numel(xp);
nm   = numel(xm);
dim1 = [nm np];
dim2 = [nm np];
A    = block_matrix(dim1,dim2);
A    = block_matrix_insert(A,dim1,dim2,1,2,Dm);
A    = block_matrix_insert(A,dim1,dim2,2,1,Dp);



if SBP_WEAK
  % SAT terms
  e0p = spalloc(np,1,1);
  e0m = spalloc(nm,1,1);
  enp = spalloc(np,1,1);
  enm = spalloc(nm,1,1);
  e0p(1)   = 1; e0m(1)   = 1;
  enp(end) = 1; enm(end) = 1;

  E01 =   -  Pmi*Z*e0m*e0m';
  En1 =   -  Pmi*Z*enm*enm';
  E02 =      Ppi*e0p*e0m';
  En2 =   -  Ppi*enp*enm';

  Dp = Dp + E02 + En2;

  if RK4
    E = block_matrix(dim1,dim2);
    E = block_matrix_insert(E,dim1,dim2,1,1,E01+En1);
    E = block_matrix_insert(E,dim1,dim2,2,1,E02+En2);
    A = A + E;
  end

end

gaussian = @(x) (x-0.5).^2.*exp(-a*(x-0.5).^2);
gaussian = @(x) exp(-a*(x-0.5).^2);

% Intialize solution
v = gaussian(xp);
if LF
  s = gaussian(xm+0.50*dt/c);
else
  s = gaussian(xm);
end


q = [s;v];
g = @(q,t) A*q;

nt = ceil(T/dt);

t = 0;


for i=1:nt

  if RK4
    q = lsrk4(g,q,t,dt);
  end

  % Leap-frog
  if LF
    v   = v + rhoi*dt*Dp*s;
    s   = s + mu*dt*Dm*v;
    if SBP_WEAK
      s = s + mu*dt*(E01 + En1)*s;
    end
  end

  t = i*dt;

  if show_plot
    plot(xm,s,'b',xp,v,'r',xm,gaussian(xm));
    ylim([-1 1]);
    pause;
  end
end  
if RK4
  s = q(1:nm); v = q(nm+1:end);
end
if ~show_plot
  plot(xm,s,'b',xm,gaussian(xm+0.5*dt/c),'k-',xp,v,'r',xp,gaussian(xp),'g');
end

% Compute error
if LF
  err.s  = norm(s  - gaussian(xm+0.5*dt/c))/norm(gaussian(xm+0.5*dt/c));
else
  err.s  = norm(s  - gaussian(xm))/norm(gaussian(xm));
end
err.v  = norm(v  - gaussian(xp))/norm(gaussian(xp));
end



