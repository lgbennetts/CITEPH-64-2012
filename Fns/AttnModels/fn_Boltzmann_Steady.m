% function fn_Boltzmann_Steady
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

function out = fn_Boltzmann_Steady(fortyp, lam0, conc, ...
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

if ~exist('DTYP','var'); DTYP=0; end
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
 cprintf([0.3,0.3,0.3],'<-------- Boltzmann steady ------->\n')
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

%% Discrete system
% Boltzmann eqn: cos(theta)*dI/dx = -beta*I + int_{-pi}^{pi} S(th,th')*I(th') dth'

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

%  N = (length(S_Fou(:,1))-1)/2;
% 
%  L_mat = zeros(2*N+1); R_mat1 = L_mat;
% 
%  L_mat(1,2)=0.5; L_mat(2*N+1,2*N)=0.5;
% 
%  for loop_N=2:2*N
%   L_mat(loop_N,loop_N-1)=0.5; L_mat(loop_N,loop_N+1)=0.5;
%  end
% 
%  R_mat0 = -beta*ones(1,2*N+1);
% 
%  count_n=-N;
%  for loop_N=1:2*N+1
%   count_p=-N;   
%   for loop_in=1:2*N+1
%    count_q = count_p-count_n;
%    if and(count_q>=-N,count_q<=N)
%     R_mat1(loop_N,loop_N) = R_mat1(loop_N,loop_N) + ...
%       S_Fou(loop_in,N+1+count_q);
%    end
%    count_p=count_p+1;
%   end
%   count_n=count_n+1;
%  end
%  
%  R_mat = conc*(diag(R_mat0)+2*pi*R_mat1)/pi/(GeomDisk(3)^2);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 2. Eigenvalues & eigenvectors %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[V,D] = eig(R_mat,L_mat); D = diag(D);

% min(abs(D))
% abs(det(V))

%plot(real(D),imag(D),'bx')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 3. Boundary conditions %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~DTYP
 %%%%%%%%%%%%%%%%%%% 
 %%%  Inf width  %%%
 %%%%%%%%%%%%%%%%%%%
 if strcmp(wth,'inf')
  
  %%%%%%%%%%%%%%%
  if absorb==0 %%
  %%%%%%%%%%%%%%%
  
  [Im,~,Iz] = fn_ArrangeEvals(D,0);
  
  if length([Im;Iz])~=length(incs)
   cprintf('red',['Check evals' '\n'])
  end
  
  % Rearrange
  V0 = V(incs,[Im;Iz]);
  
  %%%%%%%
  else %%
  %%%%%%% 
  
   Im = fn_ArrangeEvals(D,1);
  
   if length([Im;Iz])~=length(incs)
    cprintf('red',['Check evals' '\n'])
   end
  
   % Rearrange
   V0 = V(incs,Im);
  
  end % END IF absorb==0
  
  eval(['I0 = ' fn_inc '; I0(refs)=[];'])
  
  I0 = reshape(I0,length(I0),1);
  
  c0 = V0\I0;
  
  Vx = V(:,[Im;Iz]);
  
 %%%%%%%%%%%%%%%%%%%% 
 %%% Finite width %%%
 %%%%%%%%%%%%%%%%%%%%
 else 
     
  %%% Incident wave fields
  eval(['I0 = ' fn_inc '; I0(refs)=0;'])
  eval(['I1 = 0*' fn_inc '; I1(incs)=0;'])
  Ibc = I0+I1;
  Ibc = reshape(Ibc,length(Ibc),1);
  
  %%%%%%%%%%%%%%%
  if absorb==0 %%
  %%%%%%%%%%%%%%% 
   
  %%% S(x)=b0*(x*v0+u0) + c0*v0 + Sum c+(n)*exp(D+(n)*x)*v+ + Sum c-(n)*exp(D-(n)*(x-w))*v-
  %%%                              n                           n
  %%%
  %%% where D+=[Ip] and D-=[Im]
  %%%
  %%% v0 regular eigenvector corresponding to 0 eval (multiplicity 2)
  %%% u0 generalised evec:
  %%%                      (R_mat-0*L_mat)*u0=L_mat*v0
  %%%
  %%% SL(0)=S0 and SR(w)=S1
  %%%
  %%% where SL=S[incs] and SR=S[refs] 
  
  [Im,Ip,Iz0,Iz1] = fn_ArrangeEvals(D,0);
  
  n0=length([Im;Iz0]); n1=length([Ip;Iz1]);
  
  if or(n0~=length(incs),n1~=length(refs))
   cprintf('magenta',['>>> Check evals: ' fortyp '=' num2str(lam0) '\n'])
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%% Generalised evalue for lam=0 %%%
  [VR,DR]=eig(R_mat); DR=diag(DR);
  
  % Find zero eval
  [~,i0]=min(abs(DR)); 
  
  % Solve R_mat*u0=L_mat*v0 => DR*(VR\u0)=VR\L_mat*v0;
  %
  % Pick 1st zero eval of R_mat-lam*L_mat (arbitrary)
  
  v0 = V(:,Iz0);
  
  u0=0*v0;
  
  dum_v = VR\L_mat*v0;
  
  for lp=1:length(u0)
   if lp~=i0
    u0(lp)=dum_v(lp)/DR(lp);
   end
  end
  
  u0 = VR*u0;
  
  %plot(th_vec,R_mat*u0-L_mat*v0)
  
  clear dum_v VR DR
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  %%% Rearrange
  V = V(:,[Iz0,Im.',Ip.']);
  D0 = D([Im]); D1 = D([Ip]);
    
  Vbc = zeros(length(th_vec));
  %%% BCs for x=0
  EDbc = [1;1;exp(D0*0);exp(-D1*wth)];
  Vbc(incs,:) = [u0(incs),V(incs,:)]*diag(EDbc); clear EDbc
  
  %%% BCs for x=wth
  EDbc = [1;1;exp(D0*wth);exp(D1*0)];
  Vbc(refs,:) = [wth*v0(refs)+u0(refs),V(refs,:)]*diag(EDbc); clear EDbc
  
  c = Vbc\Ibc; clear Vbc Ibc
  
  %%%%%%%
  else %%
  %%%%%%% 
  
  %%% S(x)=Sum c+(n)*exp(D+(n)*x)*v+ + Sum c-(n)*exp(D-(n)*(x-w))*v-
  %%%       n                             n
  %%%
  %%% where D+=[Im] and D-=[Ip]
  %%%
  %%% SL(0)=S0 and SR(w)=S1
  %%%
  %%% where SL=S[incs] and SR=S[refs] 
  
  [Im,Ip] = fn_ArrangeEvals(D,1);
  
  n0=length([Im]); n1=length([Ip]);
  
  if or(n0~=length(incs),n1~=length(refs))
   cprintf('magenta',['>>> Check evals: ' fortyp '=' num2str(lam0) '\n'])
  end

  %%% Rearrange
  V = [V(:,[Im]),V(:,[Ip])];
  D0 = D([Im]); D1 = D([Ip]);
    
  Vbc = zeros(length(th_vec));
  %%% BCs for x=0
  EDbc = [exp(D0*0);exp(-D1*wth)];
  Vbc(incs,:) = V(incs,:)*diag(EDbc); clear EDbc
  
  %%% BCs for x=wth
  EDbc = [exp(D0*wth);exp(D1*0)];
  Vbc(refs,:) = V(refs,:)*diag(EDbc); clear EDbc
  
  c = Vbc\Ibc; clear Vbc Ibc
  
  end % END IF absorb==0
  
 end % END IF wdt==inf

else
    
 cprintf('red',['Not coded yet' '\n'])
 
end

%%%%%%%%%%%%%%
%% 4. Plot? %%
%%%%%%%%%%%%%%

if PLOT
 if ~DTYP
    
  figure(fig); hold on %h1 = subplot(1,1,1);
  %%%%%%%%%%%%%%%%%%% 
  %%%  Inf width  %%%
  %%%%%%%%%%%%%%%%%%% 
  if strcmp(wth,'inf')
   
   %%% PLOT I(theta) for increasing values of x
   if PLOT == 1
   
   %%%%%%%%%%%%%%%
   if absorb==0 %%
   %%%%%%%%%%%%%%% 
   
   for loop_x=1:length(x)
    I = Vx*diag(exp(D([Im;Iz])*x(loop_x)))*c0;
    if ~isempty(find(abs(imag(I))>tol))
     cprintf('red',['Check I: ', int2str(loop_x), ', ', ...
      num2str(max(abs(imag(I)))), '\n'])
    end
    I = [I(end);I];
    plot(gca,[-pi,th_vec]/pi,real(I))
    set(gca, 'ylim',[0,max(I0)],'xlim',[-1,1])
    title(['x=' num2str(x(loop_x))])
    xlabel('\theta','fontsize',16); ylabel('I','fontsize',16);
    if PS==0
     cprintf('m',['>> paused: hit any key to continue \n'])
     pause
    else
     pause(PS)
    end
   end
   close(gcf)
   if COMM; cprintf('blue',['>> const = ' num2str(I(1)) '\n']); end
   
   %%%%%%%
   else %%
   %%%%%%%
   
   for loop_x=1:length(x)
    I = Vx*diag(exp(D(Im)*x(loop_x)))*c0;
    if ~isempty(find(abs(imag(I))>tol))
     cprintf('red',['Check I: ', int2str(loop_x), ', ', ...
      num2str(max(abs(imag(I)))), '\n'])
    end
    I = [I(end);I];
    plot(gca,[-pi,th_vec]/pi,real(I))
    set(gca, 'ylim',[0,max(I0)],'xlim',[-1,1])
    title(['x=' num2str(x(loop_x))])
    xlabel('\theta','fontsize',16); ylabel('I','fontsize',16);
    if PS==0
     cprintf('m',['>> paused: hit any key to continue \n'])
     pause
    else
     pause(PS)
    end
   end
   close(gcf)
   if COMM; cprintf('blue',['>> const = ' num2str(I(1)) '\n']); end
   
   end % END IF absorb ==0
   
   %%%%%%%%%%%%%%%%%%%
   elseif PLOT == 2 %%
   %%%%%%%%%%%%%%%%%%% 
   
   %%%%%%%%%%%%%%%
   if absorb==0 %%
   %%%%%%%%%%%%%%% 
    I0 = Vx*diag(exp(D([Im;Iz])*0))*c0;
    H0 = real(sqrt(sum([I0(end);I0]))); clear I0
    H_vec = 0*x;
    for loop_x=1:length(x)
     I = Vx*diag(exp(D([Im;Iz])*x(loop_x)))*c0;
     I = [I(end);I];
     H_vec(loop_x) = real(sqrt(sum(I)))/H0; clear I
    end
    plot(x,H_vec,col); set(gca,'box','on')
    xlabel('x','fontsize',14); ylabel('H_{s}(x)/H_{s}(0)','fontsize',14);
   %%%%%%%
   else %%
   %%%%%%% 
    I0 = Vx*diag(exp(D([Im])*0))*c0;
    H0 = real(sqrt(sum([I0(end);I0]))); clear I0
    H_vec = 0*x;
    for loop_x=1:length(x)
     I = Vx*diag(exp(D([Iz])*x(loop_x)))*c0;
     I = [I(end);I];
     H_vec(loop_x) = real(sqrt(sum(I)))/H0; clear I
    end
    plot(x,H_vec,col); set(gca,'box','on')
    xlabel('x','fontsize',14); ylabel('H_{s}(x)/H_{s}(0)','fontsize',14);
   end % END IF absorb==0
   end % END if PLOT==
   
  %%%%%%%%%%%%%%%%%%%% 
  %%% Finite width %%%
  %%%%%%%%%%%%%%%%%%%%
  
  else
   %%% PLOT I(theta) for increasing values of x
   if PLOT==1
   %%%%%%%%%%%%%%% 
   if absorb==0 %% 
   %%%%%%%%%%%%%%%
   for loop_x=1:length(x)
    I = [x(loop_x)*v0+u0,V]*diag([1;1;exp(D0*x(loop_x));exp(D1*(x(loop_x)-wth))])*c;
    if ~isempty(find(abs(imag(I))>tol))
     cprintf('red',['Check I: ', int2str(loop_x), ', ', ...
      num2str(max(abs(imag(I)))), '\n'])
    end
    I = [I(end);I];
    plot(gca,[-pi,th_vec]/pi,real(I))
    set(gca, 'ylim',[0,max(I0)],'xlim',[-1,1])
    title(['x=' num2str(x(loop_x))])
    xlabel('\theta','fontsize',16); ylabel('I','fontsize',16);
    %set(gca,'yscale','log')
    if PS==0
     cprintf('m',['>> paused: hit any key to continue \n'])
     pause
    else
     pause(PS)
    end
   end
   %%%%%%%
   else %%
   %%%%%%% 
   for loop_x=1:length(x)
    I = V*diag([exp(D0*x(loop_x));exp(D1*(x(loop_x)-wth))])*c;
    if ~isempty(find(abs(imag(I))>tol))
     cprintf('red',['Check I: ', int2str(loop_x), ', ', ...
      num2str(max(abs(imag(I)))), '\n'])
    end
    I = [I(end);I];
    plot(gca,[-pi,th_vec]/pi,real(I))
    set(gca, 'ylim',[0,max(I0)],'xlim',[-1,1])
    title(['x=' num2str(x(loop_x))])
    xlabel('\theta','fontsize',16); ylabel('I','fontsize',16);
    %set(gca,'yscale','log')
    if PS==0
     cprintf('m',['>> paused: hit any key to continue \n'])
     pause
    else
     pause(PS)
    end
   end
   
   end % END IF absorb==0
   
   close(gcf)
   
   % Plot Hs = 4*sqrt(m0) where m0 = int I(theta) dtheta
   % Nb. actually plot ratio Hs(x)/Hs(0)
   elseif PLOT == 2
    %%%%%%%%%%%%%%%
    if absorb==0 %%
    %%%%%%%%%%%%%%%
     I0 = [u0,V]*diag([1;1;exp(D0*0);exp(D1*(0-wth))])*c;
     H0 = real(sqrt(sum([I0(end);I0]))); clear I0
     H_vec = 0*x;
     for loop_x=1:length(x)
      I = [x(loop_x)*v0+u0,V]*diag([1;1;exp(D0*x(loop_x));exp(D1*(x(loop_x)-wth))])*c;
      I = [I(end);I];
      H_vec(loop_x) = real(sqrt(sum(I)))/H0; clear I
     end
     plot(x,H_vec,col); set(gca,'box','on')
     xlabel('x','fontsize',14); ylabel('H_{s}(x)/H_{s}(0)','fontsize',14);
    %%%%%%%
    else %%
    %%%%%%%
    
    I0 = V*diag([exp(D0*0);exp(D1*(0-wth))])*c;
    H0 = real(sqrt(sum([I0(end);I0]))); clear I0
    H_vec = 0*x;
    for loop_x=1:length(x)
     I = V*diag([exp(D0*x(loop_x));exp(D1*(x(loop_x)-wth))])*c;
     I = [I(end);I];
     H_vec(loop_x) = real(sqrt(sum(I)))/H0; clear I
    end
    plot(x,H_vec,col); set(gca,'box','on')
    xlabel('x','fontsize',14); ylabel('H_{s}(x)/H_{s}(0)','fontsize',14);
    end % END IF absorb==0
   end % END IF PLOT==
  end % ENF IF wth==inf
   
 else
  cprintf('red',['Not coded yet' '\n'])
 end  
 hold off
end

%%%%%%%%%%%%%%%%
%% 5. Outputs %%
%%%%%%%%%%%%%%%%

jj0 = find(~or(th_vec(incs)>pi*th0,th_vec(incs)<-pi*th0));

%%%%%%%%%%%%%%%%%%% 
%%%  Inf width  %%%
%%%%%%%%%%%%%%%%%%% 
if strcmp(wth,'inf')
 I0 = (V(:,[Im;Iz])*c0).';
 I  = (V(:,[Im;Iz])*diag(exp(D([Im;Iz])*Param.MIZ_length))*c0).';
 if ~isempty(find(abs(imag([I,I0]))>tol))
  cprintf('red',['Check I and I0: ', int2str(loop_x), ', ', ...
   num2str(max(abs(imag(I)))), '\n'])
 end
 j0 = find(th_vec(incs)==0);
 T0 = real(I(incs(j0)))/real(I0(incs(j0)));
 Tx = sum(cos(th_vec(incs)).*real(I(incs)))/...
  sum(cos(th_vec(incs)).*real(I0(incs)));
 Tf = sum(real(I(incs)))/sum(real(I0(incs)));
 TN = sum(real(I(incs(jj0))))/sum(real(I0(incs(jj0))));
 
%%%%%%%%%%%%%%%%%%%% 
%%% Finite width %%%
%%%%%%%%%%%%%%%%%%%% 
else
 %%%%%%%%%%%%%%%
 if absorb==0 %%
 %%%%%%%%%%%%%%%
 I0 = [0*v0(incs)+u0(incs),V(incs,:)]*diag([1;1;exp(D0*0);exp(D1*(0-wth))])*c;
 I  = [wth*v0(incs)+u0(incs),V(incs,:)]*diag([1;1;exp(D0*wth);...
  exp(D1*(wth-wth))])*c;
 if ~isempty(find(abs(imag([I;I0]))>tol))
  cprintf('red',['Check I and I0: ', int2str(loop_x), ', ', ...
   num2str(max(abs(imag(I)))), '\n'])
 end
 j0 = find(th_vec(incs)==0);
 T0 = real(I(j0))/real(I0(j0));
 Tx = sum(cos(th_vec(incs)).*real(I.'))/...
  sum(cos(th_vec(incs)).*real(I0.'));
 Tf = sum(real(I.'))/sum(real(I0.'));
 TN = sum(real(I(jj0).'))/sum(real(I0(jj0).'));
 
 %%%%%%%
 else %%
 %%%%%%% 
 I0 = V(incs,:)*diag([exp(D0*0);exp(D1*(0-wth))])*c;
 I  = V(incs,:)*diag([exp(D0*Param.MIZ_length);...
  exp(D1*(Param.MIZ_length-wth))])*c;
 if ~isempty(find(abs(imag([I;I0]))>tol))
  cprintf('red',['Check I and I0: ', int2str(loop_x), ', ', ...
   num2str(max(abs(imag(I)))), '\n'])
 end
 j0 = find(th_vec(incs)==0);
 T0 = real(I(j0))/real(I0(j0));
 Tx = sum(cos(th_vec(incs)).*real(I.'))/...
  sum(cos(th_vec(incs)).*real(I0.'));
 Tf = sum(real(I.'))/sum(real(I0.'));
 TN = sum(real(I(jj0).'))/sum(real(I0(jj0).'));
 end % END IF absorb==0
end % END IF wth=inf

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SUBFUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% function [Ip,Ii,Iz] = fn_ArrangeEvals(D,tol)
%
% OUTPUT
%
% Ip = positive real evals
% Im = positive real evals
% Ii = purely imaginary evals (arranged into conj pairs)
% Iz = zero evals

function [Im,IzA,Ip,IzB] = fn_ArrangeEvals(D,tol)

% Iz = find(abs(D)<tol);
% Im = find(and(~abs(D)<tol,real(D)<-tol));
% Ip = find(and(~abs(D)<tol,real(D)>tol));
% 
% if length(Iz) ~= 2
%  cprintf('magenta',['>>> Check zero evals' '\n'])
% end
% 
% IzA = Iz(1); IzB = Iz(2);
% 
% if ~isempty(find(abs(imag(D))>tol))
%  cprintf('blue',['>>> nb. imag component to evals: ' ...
%   num2str(max(abs(imag(D)))) '\n'])
% end

D0=D;
[~,IzA] = min(abs(D0));
D0(IzA) = nan;
[~,IzB] = min(abs(D0));
D0(IzB) = nan;
Ip = find(real(D0)>-imag(D0));
Im = find(real(D0)<-imag(D0));

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

%

% function I = fn_TrapRule(xx,yy)
% 
% xtil=diff(xx);
% ytil=0.5*(yy(1:end-1)+yy(2:end));
% 
% xtil = reshape(xtil,1,length(xtil));
% ytil = reshape(ytil,1,length(ytil));
% 
% I = sum(xtil.*ytil);
% 
% return