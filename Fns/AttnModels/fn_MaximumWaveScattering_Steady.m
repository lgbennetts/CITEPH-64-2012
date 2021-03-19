% function fn_MaximumWaveScattering_Steady
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
%
% L Bennetts 2013 / Adelaide
% 
% Last updated:     March 2015

function out = fn_MaximumWaveScattering_Steady(fortyp, lam0, conc, ...
 Param, outputs, COMM, PLOT)

    
if ~exist('PLOT','var'); PLOT=1; end
if ~exist('COMM','var'); COMM=1; end
if ~exist('PS','var');   PS  =0; end

if ~exist('outputs','var'); outputs = 'transmitted energy'; end

if ~exist('th0','var'); th0=0.025; end

%% Define test

if ~exist('fortyp','var'); fortyp='freq'; end
if ~exist('lam0','var'); lam0=1/.65; end

if ~exist('RIGID','var'); RIGID=4; end
if ~exist('SURGE','var'); SURGE=1; end

if ~exist('Vert_Modes','var'); Vert_Modes=100; end

if ~exist('conc','var'); conc=0.79; end % [0.39,0.79]

if ~exist('Param','var'); Param = ParamDef_Oceanide(RIGID);
 Param = ModParam_def(Param,1,Vert_Modes,0,0,100); end

if ~exist('fn_inc','var'); fn_inc = 'kron_delta(0,th_vec)'; end
%if ~exist('fn_inc','var'); fn_inc = 'cos(th_vec).^1000'; end
%if ~exist('fn_inc','var'); fn_inc = 'series_delta(th_vec,25)'; end

if ~exist('absorb','var'); absorb=0; end

if ~exist('wth','var'); wth=Param.MIZ_length; end

if PLOT; x = linspace(0,Param.MIZ_length,6); end

%% Numerics

if ~exist('tol','var'); tol=1e-5; end
if ~exist('TTYP','var'); TTYP='asymm'; end

th_res = Param.th_res;

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


if COMM
 cprintf([0.3,0.3,0.3],'<-------- Boltzmann with only energy loss steady ------->\n')
 cprintf([0.3,0.3,0.3],['>> ' fn_inc '\n'])
 if strcmp(wth,'inf')
  cprintf([0.3,0.3,0.3],['>> semi-inf problem \n'])
 else
  cprintf([0.3,0.3,0.3],['>> width = ' num2str(wth) '\n'])
 end
 if strcmp(fortyp,'freq')
  cprintf([0.3,0.3,0.3],['>> period = ' num2str(1/lam0) '\n'])
 else
  cprintf([0.3,0.3,0.3],['>> ' fortyp ' = ' num2str(lam0) '\n'])
 end
 cprintf([0.3,0.3,0.3],['>> ' num2str(100*conc) ' concentration \n'])
 cprintf([0.3,0.3,0.3],['>> floe diameter = ' num2str(Param.floe_diam) '\n'])
 cprintf([0.3,0.3,0.3],['>> ' num2str(Param.thickness) ' thick\n'])
 cprintf([0.3,0.3,0.3],['>> rigidity = ' sprintf('%0.5g',Param.E) '\n'])
 cprintf([0.3,0.3,0.3],['>>> Vertical modes = ' int2str(Param.Ndtm) '\n'])
 cprintf([0.3,0.3,0.3],['>>> angular reslution = ' int2str(th_res) '\n'])
end

out = fn_ElasticDisk(fortyp, lam0, Param, 'Energy', th_vec, ...
 RIGID, SURGE, 0);


for loop_out=1:length(out)
 if strcmp(out(loop_out).name,'E0')
  beta=out(loop_out).value;
 elseif strcmp(out(loop_out).name,'E')
  S=out(loop_out).value;
 end
end % end loop_out

%%%%%%%%%%%%%%%%%%%%%%%%
%% 1. Discrete system %%
%%%%%%%%%%%%%%%%%%%%%%%%

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

%% Full Solution
% Boltzmann eqn - with only energy loss: cos(theta)*dI/dx = -beta/cos(theta)*I 
% This has the solution I = I_0(theta) exp(-beta/cos(theta)x )
% Where I_0(theta), is the initial incident wave information
% Since I_0(theta) = 1/2
%%%%%%%%%%%%%%%%%%%%%%%%
c1 = 0.00212;
c2 = 0.0459;
betam = conc*(c1*lam0^2 + c2*lam0^4);

Tf = 0;
TN = sqrt(exp(-(beta + betam)*wth));
Tx = sqrt(exp(-(beta + betam)*wth));
T0 = 0;


out_str = ' ''dummy'' '; out_val = ' 0 ';

if strfind(outputs,'full trans energy')
 out_str = [out_str '; ''transmitted energy full'' '];
 out_val = [out_val '; Tf'];
end

if strfind(outputs,'trans X energy')
 out_str = [out_str '; ''transmitted energy x cos(\theta)'' '];
 out_val = [out_val '; Tx'];
end

if strfind(outputs,'transmitted energy')
 out_str = [out_str '; ''transmitted energy'' '];
 out_val = [out_val '; T0'];
end

if strfind(outputs,'int trans energy')
 out_str = [out_str '; ''transmitted energy (|\theta|<' num2str(th0) '\pi )'' '];
 out_val = [out_val '; TN'];
end

eval(['out=struct( ''name'', {' out_str ...
 '}, ''value'', {' out_val '});'])
out(1)=[];

if COMM
 cprintf([0.3,0.3,0.3],'<----- End: Boltzmann steady ----->\n')
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