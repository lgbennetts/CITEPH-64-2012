% function [I,th_vec] = Main_Boltzmann(TEST, fortyp, lam0, conc, th_res, COMM)
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

function [I,th_vec] = Main_Boltzmann(TEST, fortyp, lam0, conc, th_res, COMM, PLOT)

if ~exist('TEST','var'); TEST='Oceanide'; end
if ~exist('PLOT','var'); PLOT=1; end

%% Prelims

if ~exist('tol','var'); tol=1e-5; end

if ~exist('COMM','var'); COMM=1; end
if COMM; path(path,'../../EXTRA_MATLAB_Fns'); end
if ~exist('DTYP','var'); DTYP=0; end

if ~exist('th_res','var'); th_res=5; end
th_vec = linspace(0,1,2*th_res); th_vec = unique([-th_vec,th_vec]);
th_vec(1)=[];
refs = find(or(th_vec>0.5,th_vec<-0.5));
incs = find(~or(th_vec>0.5,th_vec<-0.5));
th_vec = pi*th_vec;

%% Define test

if strcmp(TEST,'Oceanide')
    
 if ~exist('RIGID','var'); RIGID=1; end   

 if ~exist('GeomDisk','var'); GeomDisk=[0,0,0.495,33e-3]; end
 if ~exist('Param','var'); Param = ParamDef3d_Oceanide(GeomDisk); 
    Param = ModParam_def(Param,1,10,0,0); end

 if ~exist('fortyp','var'); fortyp='waveno'; end
 if ~exist('lam0','var'); lam0=2*pi./3; end

 if ~exist('bed','var'); bed=3; end
 if ~exist('conc','var'); conc=0.79; end % [0.39,0.79]
 
 %if ~exist('fn_inc','var'); fn_inc = 'kron_delta(0,th_vec)'; end
 
 if ~exist('fn_inc','var'); fn_inc = 'cos(th_vec).^100'; end
 
 dist = 5;
 
else
    
 if ~exist('GeomDisk','var'); GeomDisk=[0,0,0.495,0.1]; end
 if ~exist('Param','var'); Param = ParamDef3d(GeomDisk); 
    Param = ModParam_def(Param,1,10,0,0); end

 if ~exist('fortyp','var'); fortyp='waveno'; end
 if ~exist('lam','var'); lam0=2*pi./5; end

 if ~exist('bed','var'); bed=2; end
 if ~exist('conc','var'); conc=0.5; end
 
 if ~exist('fn_inc','var'); fn_inc = 'cos(th_vec).^2'; end

end 
 
if ~exist('absorb','var'); absorb=0; end

Forcing = Force_def(Param.g(1), bed, fortyp, lam0);

[beta, S, S_Fou] = fn_ElasticDisk(Param, GeomDisk, Forcing, bed, th_vec, RIGID, COMM);

beta = beta + absorb;

%% Discrete system
% Boltzmann eqn: cos(theta)*dI/dx = -beta*I + int_{-pi}^{pi} S(th,th')*I(th') dth'
 
if ~DTYP
    
 % Nb. This is the trap rule (because it's periodic)
 
 dx = 2*pi/length(th_vec);

 L_mat = cos(th_vec);
 R_mat0 = -beta + 0*L_mat;
 L_mat = diag(L_mat);
 
 R_mat1 = zeros(length(th_vec));
 mk = 2*th_res+1;
 for loop_th=1:2*th_res
  %T_mat(loop_th,:) = [th_vec(mk:end),th_vec(1:mk-1)];
  R_mat1(loop_th,:) = [S(mk:end),S(1:mk-1)];   
  mk=mk-1;
 end
 %T_mat(2*th_res+1,:) = th_vec;
 R_mat1(2*th_res+1,:) = S;
 mk=length(th_vec);
 for loop_th=2*th_res+2:length(th_vec)
  %T_mat(loop_th,:) = [th_vec(mk:end),th_vec(1:mk-1)];
  R_mat1(loop_th,:) = [S(mk:end),S(1:mk-1)]; 
  mk=mk-1;
 end
 clear mk
 R_mat1 = dx*R_mat1;
 
 R_mat = conc*(diag(R_mat0)+R_mat1)/pi/(GeomDisk(3)^2);
 
else % Using Fourier series (CHECK!!!!)

 N = (length(S_Fou(:,1))-1)/2;

 L_mat = zeros(2*N+1); R_mat1 = L_mat;

 L_mat(1,2)=0.5; L_mat(2*N+1,2*N)=0.5;

 for loop_N=2:2*N
  L_mat(loop_N,loop_N-1)=0.5; L_mat(loop_N,loop_N+1)=0.5;
 end

 R_mat0 = -beta*ones(1,2*N+1);

 count_n=-N;
 for loop_N=1:2*N+1
  count_p=-N;   
  for loop_in=1:2*N+1
   count_q = count_p-count_n;
   if and(count_q>=-N,count_q<=N)
    R_mat1(loop_N,loop_N) = R_mat1(loop_N,loop_N) + ...
      S_Fou(loop_in,N+1+count_q);
   end
   count_p=count_p+1;
  end
  count_n=count_n+1;
 end
 
 R_mat = conc*(diag(R_mat0)+2*pi*R_mat1)/pi/(GeomDisk(3)^2);

end

%% Eigenvalues & eigenvectors

[V,D] = eig(R_mat,L_mat); D = diag(D);

% min(abs(D))
% abs(det(V))

%plot(real(D),imag(D),'bx')

%% Boundary conditions

if ~DTYP

 [Im,Iz] = fn_ArrangeEvals(D,tol);
 
 if length([Im;Iz])~=length(incs)
  cprintf('red',['Check evals' '\n'])
 end
 
 V0 = V(incs,[Im;Iz]);
 
 eval(['I0 = ' fn_inc '; I0(refs)=[];'])
 
 I0 = reshape(I0,length(I0),1);
 
 c0 = V0\I0;

else
    
 cprintf('red',['Not coded yet' '\n'])
 
end

Vx = V(:,[Im;Iz]);
 
%% Energy conservation
   
% I = Vx*c0;
% 
% EnInc = dx*sum(I(incs));
% EnRef = dx*sum(I(refs));
% 
% I = V(1,Iz)*c0(end);
%    
% EnTra = dx*sum(I);
% 
% if abs(-EnInc+EnRef+EnTra)>tol
%  cprintf('red',['energy error: ' num2str(abs(EnIn-EnOut)) '\n'])
% end

%% Plot?

if PLOT
 if ~DTYP
    
  figure; h1 = subplot(1,1,1);
 
  x = 0:500:5000;
 
  %Vx = V(:,[Im;Iz]);
  
  for loop_x=1:length(x)
   I = Vx*diag(exp(D([Im;Iz])*x(loop_x)))*c0;
   if ~isempty(find(abs(imag(I))>tol))
    cprintf('red',['Check I: ', int2str(loop_x), ', ', ...
       num2str(max(abs(imag(I)))), '\n'])
   end
   I = [I(end);I];
   plot(h1,[-pi,th_vec]/pi,real(I))
   set(h1, 'ylim',[0,1]); set(h1,'xlim',[-1,1])
   title(['x=' num2str(x(loop_x))])
   pause
  end
  close(gcf)
  disp(['const=' num2str(I(1))])
 else
  cprintf('red',['Not coded yet' '\n'])
 end    
end

%% Outputs

if ~PLOT

 incs = find(~or(th_vec>0.5*pi,th_vec<0));

 I = V(:,[Im;Iz])*diag(exp(D([Im;Iz])*dist))*c0;
 I = I(incs); 

 if max(abs(imag(I)))>tol
  cprintf('r',['check solution: ' num2str(max(abs(imag(I)))) '\n'])    
 end

 I = real(I);

 th_vec = th_vec(incs);
 
else
  
 I=[]; th_vec=[];

end
 
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

