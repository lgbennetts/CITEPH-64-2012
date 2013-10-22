% function [I,th_vec] = fn_Boltzmann_Steady(TEST, fortyp, lam0, conc, th_res, COMM)
%
% INPUTS:
%
% fortyp = 'freq' or 'wlength' or 'waveno'
% lam0 = value of forcing fortyp
% GeomDisk = [x-location y-location Radius thickess] (coords not needed!)
% file_marker = string for file identifier; use 0 for no write
% tol = tolerance on eigenvalues being real/imag
%
% FLAGS:
%
% DTYP = flag for discretisation type: point-wise (0) or Fourier (1)
% COMM = flag for comments on (1) or off (0)
% PLOT = flag for plots on (1) or off (0)
% SURGE = include surge motion
% RIGID = rigid disk (inf) or elastic disk (rigidity=10^RIGID)

function [I,th_vec] = fn_Boltzmann_Steady(fortyp, lam0, conc, ...
 Param, th_res, COMM, PLOT)

if ~exist('PLOT','var'); PLOT=1; end

if PLOT; x = 0:50:500; end

if ~exist('COMM','var'); COMM=1; end

if ~exist('th_res','var'); th_res=100; end

%% Prelims

if ~exist('tol','var'); tol=1e-5; end

if ~exist('DTYP','var'); DTYP=0; end
if ~exist('TTYP','var'); TTYP='asymm'; end

if strcmp(TTYP,'asymm')
 th_vec = linspace(0,1,2*th_res); th_vec = unique([-th_vec,th_vec]);
 th_vec(1)=[]; 
elseif strcmp(TTYP,'symm')
 th_vec = linspace(0,1,2*th_res+1); 
 th_vec = 0.5*(th_vec(2:end)+th_vec(1:end-1));
 th_vec = [-fliplr(th_vec),th_vec];
%  th_vec = linspace(-1,1,4*th_res); 
%  th_vec = 0.5*(th_vec(2:end)+th_vec(1:end-1));
end
 
refs = find(or(th_vec>0.5,th_vec<-0.5));
incs = find(~or(th_vec>0.5,th_vec<-0.5));
th_vec = pi*th_vec;

%% Define test

if ~exist('RIGID','var'); RIGID=10; end
if ~exist('SURGE','var'); SURGE=0; end

if ~exist('Vert_Modes','var'); Vert_Modes=1; end

if ~exist('conc','var'); conc=0.79; end % [0.39,0.79]

if ~exist('Param','var'); Param = ParamDef_Oceanide(RIGID);
 Param = ModParam_def(Param,1,Vert_Modes,0,0); end

if ~exist('fortyp','var'); fortyp='freq'; end
if ~exist('lam0','var'); lam0=1/.65; end

%if ~exist('fn_inc','var'); fn_inc = 'kron_delta(0,th_vec)'; end
if ~exist('fn_inc','var'); fn_inc = 'cos(th_vec).^100'; end
%if ~exist('fn_inc','var'); fn_inc = 'series_delta(th_vec,25)'; end

if ~exist('absorb','var'); absorb=0; end

out = fn_ElasticDisk(fortyp, lam0, Param, 'Energy', th_vec, ...
 RIGID, SURGE, 0);

for loop_out=1:length(out)
 if strcmp(out(loop_out).name,'E0')
  beta=out(loop_out).value;
 elseif strcmp(out(loop_out).name,'E')
  S=out(loop_out).value;
 end
end % end loop_out
 
%% Discrete system
% Boltzmann eqn: cos(theta)*dI/dx = -beta*I + int_{-pi}^{pi} S(th,th')*I(th') dth'
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% FINITE DIFFERENCE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~DTYP
    
 % Nb. This is the trap rule (because it's periodic)
 
 dx = 2*pi/length(th_vec);

 L_mat = cos(th_vec);
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
 
 if 0
  R_mat0 = -beta + absorb + 0*S;
 else
  R_mat0 = -sum(R_mat1,1) + absorb;
 end
 R_mat = conc*(diag(R_mat0)+R_mat1)/pi/((Param.floe_diam/2)^2);
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% FOURIER SERIES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 
else 
 
 cprintf('green','!!!This code has not been checked!!!\n')

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
 
  %Vx = V(:,[Im;Iz]);
  
  for loop_x=1:length(x)
   I = Vx*diag(exp(D([Im;Iz])*x(loop_x)))*c0;
   if ~isempty(find(abs(imag(I))>tol))
    cprintf('red',['Check I: ', int2str(loop_x), ', ', ...
       num2str(max(abs(imag(I)))), '\n'])
   end
   I = [I(end);I];
   plot(h1,[-pi,th_vec]/pi,real(I))
   set(h1, 'ylim',[0,max(I0)]); set(h1,'xlim',[-1,1])
   title(['x=' num2str(x(loop_x))])
   xlabel('\theta','fontsize',16); ylabel('I','fontsize',16); 
   pause(0.5)
   end
  close(gcf)
  if COMM; disp(['const=' num2str(I(1))]); end
 else
  cprintf('red',['Not coded yet' '\n'])
 end    
end

%% Outputs

if ~PLOT

 incs = find(~or(th_vec>0.5*pi,th_vec<0));

 I = V(:,[Im;Iz])*diag(exp(D([Im;Iz])*Param.MIZ_length))*c0;
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

% sum_{n=-N}^{N} exp{i*n*theta}
% ->delta(theta-0) as N->infty

function out=series_delta(theta,N)

out=ones(1,length(theta));

for lp=1:N
    
 out = out + 2*cos(lp*theta);
 
end

return

