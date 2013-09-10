% function [I,th_vec] = fn_Boltzmann_Time(TEST, fortyp, lam0, conc, th_res, COMM)
%
% INPUTS:
%
% fortyp = 'freq' or 'wlength' or 'waveno'
% lam0 = value of forcing fortyp
% GeomDisk = [x-location y-location Radius thickess] (coords not needed!)
% file_marker = string for file identifier; use 0 for no write
% tol = tolerance on eigenvalues being real/imag
% DTYP = flag for discretisation type: point-wise (0) or Fourier (1)
% COMM = flag for comments on (1) or off (0)
% PLOT = flag for plots on (1) or off (0)

function fn_Boltzmann_Time(TEST, PLOT, th_res)

path(path,'../../EXTRA_MATLAB_Fns');

if ~exist('TEST','var'); TEST='Oceanide'; end
if ~exist('PLOT','var'); PLOT=1; end
if ~exist('PRB','var'); PRB='semiinf'; end
if ~exist('COMM','var'); COMM=1; end
if ~exist('DTYP','var'); DTYP=0; end

%% Prelims & define test

if ~exist('tol','var'); tol=1e-5; end

if strcmp(TEST,'Oceanide')
    
 if ~exist('RIGID','var'); RIGID=1; end   

 if ~exist('GeomDisk','var'); GeomDisk=[0,0,0.495,33e-3]; end
 if ~exist('Param','var'); Param = ParamDef3d_Oceanide(GeomDisk); 
    Param = ModParam_def(Param,1,10,0,0); end

 if ~exist('fortyp','var'); fortyp='waveno'; end
 if ~exist('lam0','var'); lam0=2*pi./5; end

 if ~exist('bed','var'); bed=3; end
 if ~exist('conc','var'); conc=0.79; end % [0.39,0.79]
 
 if ~exist('fn_inc','var'); fn_inc = 'fn_incA'; end
 
else
    
 cprintf('green','put something here\n')

end 
 
if ~exist('absorb','var'); absorb=0; end

Forcing = Force_def(Param.g(1), bed, fortyp, lam0);

if ~exist('th_res','var'); th_res=3; end
% th_vec = linspace(0,1,2*th_res); th_vec = unique([-th_vec,th_vec]);
% th_vec(1)=[]; th_vec = pi*th_vec;
th_vec = linspace(-1,1,4*th_res);
dth = 2*pi/length(th_vec);
th_vec = (th_vec(1:end-1) + th_vec(2:end))/2; th_vec = pi*th_vec;

if ~exist('x_vec','var'); x_vec = -9.5:1:99.5; end
dx = x_vec(2)-x_vec(1);

if ~exist('cg','var'); cg = fn_grpvel(Forcing.kappa,bed,Forcing.f); end
if ~exist('Dt','var'); Dt = 5*dx/cg; end
if ~exist('Tsteps','var'); Tsteps = 10; end

%% Kernel: single scatterer (disk)

[beta, S] = fn_ElasticDisk(Param, GeomDisk, Forcing, bed, th_vec, RIGID, COMM);

%beta = dth*sum(S(1,:));

beta = beta + absorb;

%% Discrete system
 
% Put in form: -i*(d/dt)I = i*(d/dx)*D_mat*I + i*A_mat*I
% where I is a vector discretised in angle

if ~DTYP
    
 % Nb. This is the trap rule (because it's periodic)
 
 D_mat = cos(th_vec);
 A_mat0 = -beta + 0*D_mat;
 
 A_mat1 = zeros(length(th_vec));
 mk = 2*th_res+1;
 for loop_th=1:2*th_res
  %T_mat(loop_th,:) = [th_vec(mk:end),th_vec(1:mk-1)];
  A_mat1(loop_th,:) = [S(mk:end),S(1:mk-1)];   
  mk=mk-1;
 end
 %T_mat(2*th_res+1,:) = th_vec;
 A_mat1(2*th_res+1,:) = S;
 mk=length(th_vec);
 for loop_th=2*th_res+2:length(th_vec)
  %T_mat(loop_th,:) = [th_vec(mk:end),th_vec(1:mk-1)];
  A_mat1(loop_th,:) = [S(mk:end),S(1:mk-1)]; 
  mk=mk-1;
 end
 clear mk
 A_mat1 = dth*A_mat1;
 
 A_mat = -conc*(diag(A_mat0)+A_mat1)/pi/(GeomDisk(3)^2);
 
 A_mat = cg*A_mat; D_mat = cg*D_mat;
 
else % Using Fourier series (CHECK!!!!)

 cprintf('green','Needs coding\n')

end

%% Eigenvalues & eigenvectors

[V,Lam] = eig(A_mat); Lam = diag(Lam);

%Lam = 0*Lam;    

A_op = V*diag(exp(-Lam*Dt))/V;

if max(max(abs(imag(A_op))))>tol
 cprintf('m','check A_op\n')
end

A_op=real(A_op);

if strcmp(PRB,'inf')

%% Infinite domain prb

%%% Initialise

I = zeros(length(th_vec),length(x_vec));

for loop_th=1:length(th_vec)
 xx0(loop_th,:) = x_vec+9;   
 I(loop_th,:) = feval(fn_inc,th_vec(loop_th),xx0(loop_th,:));
end

if PLOT
 if ~DTYP
  figure;  
  surf(x_vec,th_vec/pi,I)
  shading interp
  set(gca,'zlim',[0,1])
  pause
 end % end if PLOT
end % end if ~DTYP

M = eye(length(th_vec));

for loop_t=1:Tsteps
    
 M = A_op*M;   
 
 for loop_th=1:length(th_vec)   
     
  xx0(loop_th,:) = xx0(loop_th,:) - D_mat(loop_th)*Dt;    

  I(loop_th,:) = feval(fn_inc,th_vec(loop_th),xx0(loop_th,:));
  
 end % end loop_th
     
 I = M*I;   
 
 if PLOT
  if ~DTYP
   surf(x_vec,th_vec/pi,I)
   shading interp
   set(gca,'zlim',[0,1])
   pause
  end % end if PLOT
 end % end if ~DTYP
 
end % end loop_t

if PLOT; close(gcf); end

elseif strcmp(PRB,'semiinf')
    
%% Semi-inf domain prb (x>0 scattering medium)

x=sym('x'); th=sym('th');

%I0 = symfun(fn_incA(x,th),[x th]);

I0 = symfun(exp(-10*(th^2))*exp(-((x-0).^2)/10),[x th]);
for loop_th=1:length(th_vec)
 I_vec{loop_th} = symfun(I0(x,th_vec(loop_th)),x);
end

%%% Initialise

if PLOT
 I = zeros(length(th_vec),length(x_vec));
 for loop_th=1:length(th_vec)
  I(loop_th,:) = eval(I_vec{loop_th}(x_vec));
 end
end

if PLOT   
 if ~DTYP
  figure;  
  surf(x_vec,th_vec/pi,I)
  shading interp
  set(gca,'zlim',[0,1])
  title('t=0')
  view([0,90])
  pause
 end % end if PLOT
end % end if ~DTYP

for loop_t=1:Tsteps
    
 if COMM; disp(['t=' num2str(loop_t*Dt)]); end   
    
 % Step 1: advect
 
 for loop_th=1:length(th_vec)
  I_int{loop_th} = symfun(I_vec{loop_th}(x-Dt*D_mat(loop_th)),x);
 end
 
 % Step 2: scatter
 
 for loop_th=1:length(th_vec)
  fn_dum = symfun(0*x,x);   
  for loop_th2=1:length(th_vec)
   fn_dum = symfun(fn_dum(x)+...
       heaviside(x)*A_op(loop_th,loop_th2)*I_vec{loop_th2}(x),x);
  end % end loop_th2
  I_vec{loop_th} = symfun(heaviside(-x)*I_int{loop_th}(x) + ...
      fn_dum(x), x);
  clear fn_dum
 end
 
 clear I_int
 
%  for loop_th=1:length(th_vec)
%   wts_mat = fn_interpx(x_vec,Dt*D_mat(loop_th)); 
%   I_int(loop_th,:) = wts_mat*transpose(I(loop_th,:));  
%  end

 if PLOT
  for loop_th=1:length(th_vec)
   I(loop_th,:) = eval(I_vec{loop_th}(x_vec));
  end
 end

 if PLOT
  if ~DTYP
   surf(x_vec,th_vec/pi,I)
   shading interp
   set(gca,'zlim',[0,1])
   title(['t=' num2str(loop_t*Dt)])
   view([0,90])
   pause
  end % end if PLOT
 end % end if ~DTYP
 
end % end loop_t

if PLOT; close(gcf); end

end % end PRB 

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SUBFUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% function [Ip,Ii,Iz] = fn_ArrangeEvals(D,tol)
%
% OUTPUT
%
% Ip = positive real evals
% Ii = purely imaginary evals (arranged into conj pairs)
% Iz = zero evals

function [Im,Iz] = fn_ArrangeEvals(D,tol)

Iz = find(abs(D)<tol);
Im = find(and(~abs(D)<tol,real(D)<-tol));

if length(Iz) ~= 2
 cprintf('red',['Check zero evals' '\n'])
end

Iz = Iz(1);

if ~isempty(find(abs(imag(D))>tol))
 cprintf('blue',['nb. imag component to evals: ' num2str(max(abs(imag(D)))) '\n'])
end

return

% Kronecker Delta

function out=kron_delta(u,v)

out=zeros(1,length(v));

for loop=1:length(v)
    
 if u==v(loop)
  out(loop)=1;
 end
 
end

return

% cosine incident wave

function out = fn_incA(th,x)

for loop_th=1:length(th)

 if ~or(th(loop_th)>0.5*pi,th(loop_th)<-0.5*pi)

  out(loop_th,:) = (cos(th(loop_th))^2)*exp(-(x.^2)/10);
 
 else
    
  out(loop_th,:) = 0*x;
 
 end % loop_th
end % if loop_th...

return

% Group velocity

function cg = fn_grpvel(k,h,f)

% see Mei book p.18

cg = 0.5*(2*pi*f/k)*(1+(2*k*h/sinh(2*k*h)));

return

% 

function wts = fn_interpx(xx,Dt)

wts = zeros(length(xx));

% assume |Dt|<=dx

if abs(Dt)>xx(2)-xx(1)
 cprintf('red','reduce time step!')   
end

ratio = abs(Dt)/(xx(2)-xx(1));

if Dt>0
 wts(1,1) = 1-ratio;
 for loop_x=2:length(xx)
   wts(loop_x,loop_x-1) = ratio; 
   wts(loop_x,loop_x) = 1-ratio;  
 end
else
 for loop_x=1:length(xx)-1
   wts(loop_x,loop_x+1) = ratio; 
   wts(loop_x,loop_x) = 1-ratio;  
 end   
 wts(length(xx),length(xx)) = 1-ratio;
end

return