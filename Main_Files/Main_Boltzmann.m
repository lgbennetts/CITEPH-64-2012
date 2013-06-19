% function Main_Boltzmann
%
% INPUTS:
%
% fortyp = 'freq' or 'wlength' or 'waveno'
% lam0 = value of forcing fortyp
% file_marker = string for file identifier; use 0 for no write
% tol = tolerance on eigenvalues being real/imag
% DTYP = flag for discretisation type: point-wise (0) or Fourier (1)
% COMM = flag for comments on (1) or off (0)

function Main_Boltzmann

%% Prelims

if ~exist('tol','var'); tol=1e-3; end

if ~exist('COMM','var'); COMM=1; end
if COMM; path(path,'../../EXTRA_MATLAB_Fns'); end
if ~exist('DTYP','var'); DTYP=0; end

if ~exist('th_res','var'); th_res=25; end
th_vec = pi*linspace(-1,1,4*th_res+1);
th_vec(4*th_res+1)=[];

if ~exist('GeomDisk','var'); GeomDisk=[0,0,0.5,0.1]; end
if ~exist('Param','var'); Param = ParamDef3d(GeomDisk); 
    Param = ModParam_def(Param,0,0); end

if ~exist('fortyp','var'); fortyp='waveno'; end
if ~exist('lam','var'); lam0=2*pi./5; end

if ~exist('bed','var'); bed=2; end
if ~exist('conc','var'); conc=0.5; end

if ~exist('absorb','var'); absorb=0; end

Forcing = Force_def(Param.g(1), bed, fortyp, lam0);

if ~exist('fn_inc','var'); fn_inc = 'cos(th_vec).^2'; end

[beta, S, S_Fou] = fn_ElasticDisk(Param, GeomDisk, Forcing, bed, th_vec, COMM);

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
 
 else

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

end

R_mat = conc*(diag(R_mat0)+2*pi*R_mat1)/pi/(GeomDisk(3)^2);

%% Eigenvalues & eigenvectors

[V,D] = eig(L_mat,R_mat); D = diag(D);

% min(abs(D))
% abs(det(V))

%plot(real(D),imag(D),'bx')

%% Boundary conditions

if ~DTYP

 eval(['I0 = ' fn_inc '; I0([1:th_res,3*th_res+2:4*th_res])=[];'])
 
 I0 = reshape(I0,2*th_res+1,1);
 
 [Im,Ii,Iz] = fn_ArrangeEvals(D,tol);
 
 inds0 = th_res+1:3*th_res+1;
 
 V0 = [V(inds0,[Im;Iz]),real(V(inds0,Ii))];
 
 c0 = V0\I0;

else
    
 cprintf('red',['Not coded yet' '\n'])
 
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

function [Im,Ii,Iz] = fn_ArrangeEvals(D,tol)

Iz = find(abs(D)<tol);
Im = find(real(D)<=-tol);
Ii(:,1) = find(imag(D)>=tol);
%Ii(:,2) = find(imag(D)<=-tol);

if length(Iz) ~= 2
 cprintf('red',['Check zero evals' '\n'])
end

if length(Ii(:,1)) ~= 2
 cprintf('red',['Check imag evals' '\n'])
end


return


