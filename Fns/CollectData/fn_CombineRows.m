% function fn_CombineRows(file_marker,Rows,ens0,row_sep,COMM)
%
% Calculates reflection/transmission by a single realisation of scatterers
% enclosed by walls.
% Location of scatterers is deterministic.
% Can be arranged into rows for numerical efficiency.
%
% INPUTS:
%
% Rows = number of rows
% ens0 = the size of the ensemble
% file_marker = string for file identifier
% RWSP = mean separation of the rows (centre-to-centre)
%           0 reverts to default perturbed periodic arrangement
% COMM = flag for comments on (1) or off (0)
%
% OUTPUTS:
%
% TraStat = array size 8xM where M is num tank modes
%           1st row contains mean energy for each mode
%           2nd row contains corresponding std
%           3rd row is energy from periodic problem
%           4th row is energy from periodic problem zeroth order
%           5th row contains mean energy for each mode zeroth order apx
%           6th row contains corresponding std
%           7th row is v_vec = cosine of angles wrt x-axis
%           8th row q's: exp(1i*q) is eval of periodic prb
%           9th corresponding zeroth order eval
% EnTStat = size 6x1 array: similar to TraStat but for energy flux x-dir 
% attn    = 6x1 array of attenuation coefficients
%           1st entry is mean attenuation rate over perturbed simulations
%              given by -log(mean trans energy)/mean array length
%           2nd entry is corresponding value for periodic problem
%           3rd entry is max_j{imag(log(D_j)/1i))}/row length
%              where D_j are evals of transfer matrix
%           4th-5th rows as in 2nd-3rd but zeroth order approx (ZEROTH=1)
%           6th entry is from problem considering interactions of 0th order
%              waves only (ZEROTH=1)

function [TraStat,EnTStat,attn]=fn_CombineRows(file_marker,Rows,ens0,RWSP,COMM)

path(path,'../../EXTRA_MATLAB_Fns')
 
%% Inputs & prelims

if ~exist('ZEROTH','var'); ZEROTH=1; end

if ~exist('ens0','var'); ens0=100; end

if ~exist('file_marker','var')
 cprintf('blue','NEED A FILE TO LOAD\n');
 file_marker = input('Which problem? ', 's');
end

if ~isempty(file_marker)
 if ~exist('COMM','var'); COMM=1; end
 
 file_pre = '../../../../../Documents/MatLab/Data/Wavetank/Rows/';
 file_name = ['Main_Rows_',file_marker,'.mat'];
 load([file_pre,file_name],'Rm_sv', 'Tm_sv', 'Rp_sv', 'Tp_sv',...
    'Rm_sv0','Tm_sv0','Rp_sv0','Tp_sv0','GeomDisks_sv0', ...
    'lam_vec','k0_sv','v_sv','x_lim_sv','epsilon','w','GeomDisks_sv',...
    'clockout','duration','v_info');

 ens = length(Rm_sv);
 Np = length(GeomDisks_sv{1}(:,1));   
 if ~exist('Rows','var'); Rows=5; end
 if ~exist('row_sep','var'); row_sep=0; end
 
%% Info?

 if COMM
  if exist('clockout','var'); 
   cprintf([0.3,0.3,0.3],['Created: ' int2str(clockout(3)),'/', int2str(clockout(2)),'/',int2str(clockout(1)),' ' ,...
    int2str(clockout(4)),':', int2str(clockout(5)),':',int2str(round(clockout(6))),'\n'] ); 
  end
  if exist('duration','var'); 
   cprintf([0.3,0.3,0.3],['Lasted: ' int2str(round(duration(1))),' hours ', ...
     int2str(round(duration(2))),' mins ',int2str(round(duration(3))) ' secs\n'] ); 
  end
  if exist('v_info','var')
   for loop=1:length(v_info)  
    cprintf([0.3,0.3,0.3],[v_info{loop,:},'\n']); 
   end
  end
 end

 %% Begin calculations: 
 
 if COMM
  cprintf([0.3,0.3,0.3],'-----------------------------------------------\n')
  cprintf([0.3,0.3,0.3],'----------   START: Combine_Rows    -----------\n')
 end
 
 for loop_lam=1:length(lam_vec)
  dum(loop_lam) = length(v_sv{loop_lam});
 end
 
 RefStat = zeros(2,max(dum),length(lam_vec)); 
 TraStat = zeros(9,max(dum),length(lam_vec));
 EnRStat = zeros(2,length(lam_vec));
 EnTStat = zeros(5,length(lam_vec));
 attn = zeros(5,length(lam_vec));
 
 clear dum
 
 %%%%%%%%%%%%%%%%%%%%
 %%% PERIODIC PRB %%%
 %%%%%%%%%%%%%%%%%%%%
 
 if ~RWSP; row_sep=w/Np; end
 sep0 = row_sep - 0.99;
   
 for loop_lam=1:length(lam_vec)

   k0 = k0_sv{loop_lam};
   v_vec = v_sv{loop_lam};
   v1 = 1:length(v_vec); v2 = v1 + length(v_vec);
   
   real_inds = find(real(v_vec)~=0); 
   epsm = ones(1,length(real_inds)); epsm(1) = 2; epsm = epsm/2;
   
   kv = reshape(diag(k0)*v_vec,1,size(v_vec,1)*size(v_vec,2));
 
   S = [Rm_sv0{1},Tp_sv0{1};Tm_sv0{1},Rp_sv0{1}];
 
   P = fn_TransScat(S,-1); 
   Ep = exp(1i*kv*sep0/2); Em = exp(-1i*kv*sep0/2);
   E = diag([Ep,Em]);
 
   [X,D] = eig(E*P*E); 
   
   [X,D,q]=fn_ArrangeTrans(X,D);
   
   invX = eye(2*length(v_vec))/X;
   
   X11=invX(v1,v1); X12=invX(v1,v2); X21=invX(v2,v1); X22=invX(v2,v2); 
   
   qN = diag(exp(1i*q(v1).*Rows));
   
   %%% Scattering matrix for multiple rows
   
   SN = -[X22, -qN*X21; qN*X12, -X11]\[X21, -qN*X22; qN*X11, -X12];
   
   TN = SN(v2,v1); T_vec = TN(real_inds,real_inds(1));
   
   TraStat(3,:,loop_lam) = abs(epsm.*transpose(T_vec)).^2;
   
   TraStat(8,:,loop_lam) = q(v1);
   
   EnTStat(3,loop_lam) = sum( transpose(epsm.*v_vec(real_inds)) .* ...
     ( abs(T_vec).^2 ) );
 
   attn(2,loop_lam) = -log(EnTStat(3,loop_lam))/Rows/row_sep;
  
   %%%% Rough estimate of attenuation rate
   attn(3,loop_lam) = max(imag(log(diag(D))/1i))/row_sep;
   
   if ZEROTH
    S = [Rm_sv0{1}(real_inds(1),real_inds(1)),...
         Tp_sv0{1}(real_inds(1),real_inds(1));...
         Tm_sv0{1}(real_inds(1),real_inds(1)),...
         Rp_sv0{1}(real_inds(1),real_inds(1))];
 
    P = fn_TransScat(S,-1); 
    Ep = exp(1i*kv(real_inds(1))*sep0/2); Em = exp(-1i*kv(real_inds(1))*sep0/2);
    E = diag([Ep,Em]);
 
    [X,D] = eig(E*P*E); 
   
    [X,D,q]=fn_ArrangeTrans(X,D);
   
    invX = eye(2)/X;
   
    v1=1; v2=2;
   
    X11=invX(v1,v1); X12=invX(v1,v2); X21=invX(v2,v1); X22=invX(v2,v2); 
   
    qN = diag(exp(1i*q(v1).*Rows));
   
    %%% Scattering matrix for multiple rows
   
    SN = -[X22, -qN*X21; qN*X12, -X11]\[X21, -qN*X22; qN*X11, -X12];
   
    T_vec = SN(v2,v1); 
     
    TraStat(4,1,loop_lam) = abs(T_vec).^2;
    
    TraStat(9,1,loop_lam) = q(1);
   
    EnTStat(4,loop_lam) = abs(T_vec).^2;
   
    attn(4,loop_lam) = -log(EnTStat(4,loop_lam))/Rows/row_sep;
  
    %%%% Rough estimate of attenuation rate
    attn(5,loop_lam) = max(imag(log(diag(D))/1i))/row_sep;
 
   end % end ZEROTH
 end
 
 clear sep0 S P Ep Em E X D q qN TN T_vec epsm real_inds
 
 %%%%%%%%%%%%%%%%%%%%
 %%% COMBINE ROWS %%%
 %%%%%%%%%%%%%%%%%%%%
 
 dist = zeros(1,ens0);

 for loop_ens=1:ens0 
  
  for loop=1:Rows; Rowperm(loop)=randi(ens); end   
     
  %%% Calculate separation of closest points of adjacent rows %%%
  %%% Assumes mean axis of rows is x=row_sep/2 %%%
  dum_dist = x_lim_sv{Rowperm(1)}(2)-x_lim_sv{Rowperm(1)}(1);
  if ~RWSP
   row_seps = zeros(1,Rows-1);
    
   for loop=1:Rows-1
    x0 = x_lim_sv{Rowperm(loop)}(2);
    x1 = x_lim_sv{Rowperm(loop+1)}(1);
   
    dum_dist = dum_dist + ...
       x_lim_sv{Rowperm(loop+1)}(2)-x_lim_sv{Rowperm(loop+1)}(1);
  
    row_seps(loop) = row_sep - (x0 - (row_sep/2)) + (x1 - (row_sep/2));
   end
   dist(loop_ens) = dum_dist + sum(row_seps); clear dum_dist
  else
   row_seps = RWSP*rand(1,Rows-1);
   % Keep distance consistent with periodic case
   dist = Rows*w/Np; 
  end
  
  fd = find(row_seps<0); count=1;
  
  while and(~isempty(fd),count<=100)
   row_seps(fd) = row_sep + 2*(rand(1,length(fd))-0.5);
   fd = find(row_seps<0);
   count = count+1;
  end
  clear count fd
  
  for loop_lam=1:length(lam_vec)
 
   %%% Frequency dependent parameters %%%
   k0 = k0_sv{loop_lam};
   v_vec = v_sv{loop_lam};
 
   %%% Initialise Scattering Matrix %%%
   %%% First Row %%%
   Rm = Rm_sv{Rowperm(1)}{loop_lam};
   Tm = Tm_sv{Rowperm(1)}{loop_lam};
   Rp = Rp_sv{Rowperm(1)}{loop_lam};
   Tp = Tp_sv{Rowperm(1)}{loop_lam};
   
   if ZEROTH
    Rm_z = Rm_sv{Rowperm(1)}{loop_lam}(1,1);
    Tm_z = Tm_sv{Rowperm(1)}{loop_lam}(1,1);
    Rp_z = Rp_sv{Rowperm(1)}{loop_lam}(1,1);
    Tp_z = Tp_sv{Rowperm(1)}{loop_lam}(1,1);
   end
   %---------------% 
   
   kv = reshape(diag(k0)*v_vec,1,size(v_vec,1)*size(v_vec,2));
 
   %%% Nb. see Main_Channel_freq for shift of array into tank %%%  
   
   v1 = 1:size(v_vec,1)*size(v_vec,2); v2 = v1+size(v_vec,1)*size(v_vec,2);
 
   Full_S = zeros(2*length(v_vec),2*length(v_vec),Rows);
 
   Full_S(v1,v1,1)=Rm; Full_S(v2,v2,1)=Rp;
   Full_S(v1,v2,1)=Tp; Full_S(v2,v1,1)=Tm;
  
   if ZEROTH
    Full_S_z = zeros(2,2,Rows);
    Full_S_z(1,1,1)=Rm_z; Full_S_z(2,2,1)=Rp_z;
    Full_S_z(1,2,1)=Tp_z; Full_S_z(2,1,1)=Tm_z;
   end
  
   %%% Iterate Rows %%%
   II = eye(size(v_vec,1)*size(v_vec,2));
   for loop_Sc=2:Rows
    % - left-to-right method
    Ep = diag(exp(1i*kv*row_seps(loop_Sc-1)));
    
    %%% Next Row %%%
    Rm = Rm_sv{Rowperm(loop_Sc)}{loop_lam};
    Tm = Tm_sv{Rowperm(loop_Sc)}{loop_lam};
    Rp = Rp_sv{Rowperm(loop_Sc)}{loop_lam};
    Tp = Tp_sv{Rowperm(loop_Sc)}{loop_lam};
    %---------------% 
    
    R0m = Full_S(v1,v1,loop_Sc-1); T0p = Full_S(v1,v2,loop_Sc-1);
    T0m = Full_S(v2,v1,loop_Sc-1); R0p = Full_S(v2,v2,loop_Sc-1);

    IV1 = II/(II - Rm*Ep*R0p*Ep); IV2 = II/(II - R0p*Ep*Rm*Ep);

    Full_S(v1,v1,loop_Sc) = R0m + T0p*Ep*IV1*Rm*Ep*T0m;
    Full_S(v1,v2,loop_Sc) = T0p*Ep*IV1*Tp;
    Full_S(v2,v2,loop_Sc) = Rp + Tm*Ep*IV2*R0p*Ep*Tp;
    Full_S(v2,v1,loop_Sc) = Tm*Ep*IV2*T0m;
    
    if ZEROTH
        
     % - left-to-right method
     Ep_z = diag(exp(1i*kv(1)*row_seps(loop_Sc-1)));
    
     %%% Next Row %%%
     Rm_z = Rm_sv{Rowperm(loop_Sc)}{loop_lam}(1,1);
     Tm_z = Tm_sv{Rowperm(loop_Sc)}{loop_lam}(1,1);
     Rp_z = Rp_sv{Rowperm(loop_Sc)}{loop_lam}(1,1);
     Tp_z = Tp_sv{Rowperm(loop_Sc)}{loop_lam}(1,1);
     %---------------% 
    
     R0m_z = Full_S_z(1,1,loop_Sc-1); T0p_z = Full_S_z(1,2,loop_Sc-1);
     T0m_z = Full_S_z(2,1,loop_Sc-1); R0p_z = Full_S_z(2,2,loop_Sc-1);

     IV1 = 1/(1 - Rm_z*Ep_z*R0p_z*Ep_z); IV2 = 1/(1 - R0p_z*Ep_z*Rm_z*Ep_z);

     Full_S_z(1,1,loop_Sc) = R0m_z + T0p_z*Ep_z*IV1*Rm_z*Ep_z*T0m_z;
     Full_S_z(1,2,loop_Sc) = T0p_z*Ep_z*IV1*Tp_z;
     Full_S_z(2,2,loop_Sc) = Rp_z + Tm_z*Ep_z*IV2*R0p_z*Ep_z*Tp_z;
     Full_S_z(2,1,loop_Sc) = Tm_z*Ep_z*IV2*T0m_z;
    end % end ZEROTH
   end % end loop_Sc
   
   %%% Collect reflection & transmission matrices %%%
   R_vec{loop_ens,loop_lam} = Full_S(v1,v1,Rows); 
   T_vec{loop_ens,loop_lam} = Full_S(v2,v1,Rows);
   clear Full_S R0m R0p T0m T0p Rm Rp Tm Tp
   if ZEROTH
    R_vec_z{loop_ens,loop_lam} = Full_S_z(1,1,Rows); 
    T_vec_z{loop_ens,loop_lam} = Full_S_z(2,1,Rows);
    clear Full_S_z R0m_z R0p_z T0m_z T0p_z Rm_z Rp_z Tm_z Tp_z
   end
   
  end % end loop_lam
 end % end loop_ens
 
 %% Outputs:
 
 %%% CHECK ENERGY ERROR %%%
 
 for loop_ens=1:ens0
  for loop_lam=1:length(lam_vec)
   v_vec = v_sv{loop_lam};
 
   disp_error = []; mkr = 1;
   real_inds = find(real(v_vec)~=0); Y0_Dim = length(real_inds);
   Vert_Dim = size(R_vec{loop_lam},1)/length(v_vec);
   epsm = ones(1,length(real_inds)); epsm(1) = 2; epsm = epsm/2;
   for loop_Y=1:Vert_Dim:Y0_Dim*Vert_Dim
    EngErr = epsm(loop_Y)*v_vec(real_inds(loop_Y)) ...
      - sum(epsm.*v_vec(real_inds).*...
       ( abs(R_vec{loop_ens,loop_lam}(real_inds,real_inds(loop_Y)).').^2 ...
       + abs(T_vec{loop_ens,loop_lam}(real_inds,real_inds(loop_Y)).').^2 ) );

    if abs(EngErr)>1e-5
     if mkr 
      disp_error{1} = ['Energy error ' 'freq = ' num2str(lam_vec(loop_lam))]; 
      disp_error{2} = '; Modes: '; disp_error{3} = ' -> '; mkr=0;
     end
     disp_error{2} = [disp_error{2} int2str(loop_Y) ', '];
     disp_error{3} = [disp_error{3} num2str(abs(EngErr))  ', '];
    end   
   end % end loop_Y
  end % end loop_lam
 end % end loop_ens
 
 if ~isempty(disp_error); cprintf('magenta', ...
  [disp_error{1},disp_error{2},disp_error{3},'\n']); end

%  for loop_ens=1:ens0
%   for loop_lam=1:length(lam_vec)
%    v_vec = v_sv{loop_lam};
%  
%    %disp_error = []; mkr = 1;
%    real_inds = find(real(v_vec)~=0); Y0_Dim = length(real_inds);
%    Vert_Dim = size(R_vec{loop_lam},1)/length(v_vec);
%    
%    for loop_Y=1:Vert_Dim:Y0_Dim*Vert_Dim
%    
%     En = sum(...
%          abs(R_vec{loop_ens,loop_lam}(real_inds,real_inds(loop_Y)).').^2 ...
%        + abs(T_vec{loop_ens,loop_lam}(real_inds,real_inds(loop_Y)).').^2 );
%    
%    end % end loop_Y
%   end % loop_lam
%  end % end loop_ens
   
 %%% DIRECTIONAL SPECTRUM REF/TRANS ENERGY %%%
 %%% NB. SOLUTION NOT NECESSARILY SYMMETRIC (PHASES MISSING HERE) %%%
 
 for loop_lam=1:length(lam_vec)
  for loop_v=1:length(v_sv{loop_lam})
      dumR=zeros(1,ens0); dumT=dumR;
   for loop_ens=1:ens0
    dumR(loop_ens) = (abs(R_vec{loop_ens,loop_lam}(loop_v,1)).^2);
    dumT(loop_ens) = (abs(T_vec{loop_ens,loop_lam}(loop_v,1)).^2);
   end
   if loop_v==1
    RefStat(1,loop_v,loop_lam) = mean(dumR);
    RefStat(2,loop_v,loop_lam) = std(dumR);
    TraStat(1,loop_v,loop_lam) = mean(dumT);
    TraStat(2,loop_v,loop_lam) = std(dumT);
   else
    RefStat(1,loop_v,loop_lam) = 0.25*mean(dumR);
    RefStat(2,loop_v,loop_lam) = 0.25*std(dumR);
    TraStat(1,loop_v,loop_lam) = 0.25*mean(dumT);
    TraStat(2,loop_v,loop_lam) = 0.25*std(dumT);
   end
  end
  TraStat(7,:) = v_vec;
 end
 
 %%% FULL ENERGY %%%
 
 for loop_lam=1:length(lam_vec)
  v_vec = v_sv{loop_lam};
  real_inds = find(real(v_vec)~=0); Y0_Dim = length(real_inds);
  epsm = ones(1,length(real_inds)); epsm(1) = 2; epsm = epsm/2;
  EnR = zeros(1,ens0); EnT=EnR;
  for loop_ens=1:ens0
   
   EnR(loop_ens) = sum( transpose(epsm.*v_vec(real_inds)) .* ...
     ( squeeze(abs(R_vec{loop_ens,loop_lam}(real_inds,real_inds(1))).^2) ) );
   EnT(loop_ens) = sum( transpose(epsm.*v_vec(real_inds)) .* ...
     ( squeeze(abs(T_vec{loop_ens,loop_lam}(real_inds,real_inds(1))).^2) ) );
 
   if ZEROTH
    EnR_z(loop_ens) = v_vec(1)*squeeze(abs(R_vec_z{loop_ens,loop_lam})).^2;
    EnT_z(loop_ens) = v_vec(1)*squeeze(abs(T_vec_z{loop_ens,loop_lam})).^2;
   end
 
  end % end loop_ens
  
  EnRStat(1,loop_lam) = mean(EnR);
  EnRStat(2,loop_lam) = std(EnR);
  EnTStat(1,loop_lam) = mean(EnT);
  EnTStat(2,loop_lam) = std(EnT);
  
  attn(1,loop_lam) = -log(EnTStat(1,loop_lam))/mean(dist);
  
  if ZEROTH
   EnTStat(5,loop_lam) = mean(EnT_z); TraStat(5,1,loop_lam) = mean(EnT_z);
   EnTStat(6,loop_lam) = std(EnT_z);  TraStat(6,1,loop_lam) = std(EnT_z);
   attn(6,loop_lam) = -log(mean(EnT_z))/mean(dist);
  end
  
 end % end loop_lam
  
 if COMM
 cprintf([0.3,0.3,0.3],'-----------    END: Combine_Rows   ------------\n')
 cprintf([0.3,0.3,0.3],'-----------------------------------------------\n')
 end
end % end if

return