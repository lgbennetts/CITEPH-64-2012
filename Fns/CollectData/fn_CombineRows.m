% function fn_CombineRows(Rows,row_sep,ens0,file_marker,COMM)
%
% Calculates reflection/transmission by a single realisation of scatterers
% enclosed by walls.
% Location of scatterers is deterministic.
% Can be arranged into rows for numerical efficiency.
%
% INPUTS:
%
% Rows = number of rows
% row_sep = mean separation of the rows (centre-to-centre)
% ens0 = the size of the ensemble
% file_marker = string for file identifier
% COMM = flag for comments on (1) or off (0)

function [TraStat,EnTStat]=fn_CombineRows(file_marker,Rows,row_sep,ens0,COMM)

if ~exist('ens0','var'); ens0=100; end

if ~exist('file_marker','var')
 path(path,'../../EXTRA_MATLAB_Fns')
 cprintf('blue','NEED A FILE TO LOAD\n');
 file_marker = input('Which problem? ', 's');
end

if ~isempty(file_marker)
 if ~exist('COMM','var'); COMM=1; end
 
 file_pre = '../../../../../Documents/MatLab/Data/Wavetank/Rows/';
 file_name = ['Main_Rows_',file_marker,'.mat'];
 load([file_pre,file_name],'Rm_sv', 'Tm_sv', 'Rp_sv', 'Tp_sv',...
    'lam_vec','k0_sv','v_sv','x_lim_sv','epsilon','w','GeomDisks_sv',...
    'clockout','duration','v_info');

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

 ens = length(Rm_sv);
 Np = length(GeomDisks_sv{1}(:,1));   
 if ~exist('Rows','var'); Rows=5; end
 if ~exist('row_sep','var'); row_sep=w/Np; end
 
 if COMM
  cprintf([0.3,0.3,0.3],'-----------------------------------------------\n')
  cprintf([0.3,0.3,0.3],'----------   START: Combine_Rows    -----------\n')
 end
 
 
 %%%%%%%%%%%%%%%%%%%%
 %%% COMBINE ROWS %%%
 %%%%%%%%%%%%%%%%%%%%

 for loop_ens=1:ens0 
  
  for loop=1:Rows; Rowperm(loop)=randi(ens); end   
     
  %%% Calculate separation of closest points of adjacent rows %%%
  %%% Assumes mean axis of rows is x=row_sep/2 %%%
  row_seps = zeros(1,Rows-1);
  for loop=1:Rows-1
   x0 = x_lim_sv{Rowperm(loop)}(2);
   x1 = x_lim_sv{Rowperm(loop+1)}(1);
   
   row_seps(loop) = row_sep - (x0 - (row_sep/2)) + (x1 - (row_sep/2));
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
   %---------------% 
   
   kv = reshape(diag(k0)*v_vec,1,size(v_vec,1)*size(v_vec,2));
 
   %%% Nb. see Main_Channel_freq for shift of array into tank %%%  
   
   v1 = 1:size(v_vec,1)*size(v_vec,2); v2 = v1+size(v_vec,1)*size(v_vec,2);
 
   Full_S = zeros(2*length(v_vec),2*length(v_vec),Rows);
 
   Full_S(v1,v1,1)=Rm; Full_S(v2,v2,1)=Rp;
   Full_S(v1,v2,1)=Tp; Full_S(v2,v1,1)=Tm;
   
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
    
   end % end loop_Sc
   
   %%% Collect reflection & transmission matrices %%%
   R_vec{loop_ens,loop_lam} = Full_S(v1,v1,Rows); 
   T_vec{loop_ens,loop_lam} = Full_S(v2,v1,Rows);
   
   clear Full_S R0m R0p T0m T0p Rm Rp Tm Tp

  end % end loop_lam
 end % end loop_ens
 
 %%%%%%%%%%%%%%
 %%% OUTPUT %%%
 %%%%%%%%%%%%%%
 
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
      disp_error{1} = ['Energy error ' 'wk/pi = ' num2str(w*lam_vec(loop_lam)/pi)]; 
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

 %%% DIRECTIONAL SPECTRUM REF/TRANS ENERGY %%%
 
 for loop_lam=1:length(lam_vec)
  for loop_v=1:length(v_sv{loop_lam})
      dumR=zeros(1,ens0); dumT=dumR;
   for loop_ens=1:ens0
    dumR(loop_ens) = abs(R_vec{loop_ens,loop_lam}(loop_v,1)).^2;
    dumT(loop_ens) = abs(T_vec{loop_ens,loop_lam}(loop_v,1)).^2;
   end
   RefStat(1,loop_v,loop_lam) = mean(dumR);
   RefStat(2,loop_v,loop_lam) = std(dumR);
   TraStat(1,loop_v,loop_lam) = mean(dumT);
   TraStat(2,loop_v,loop_lam) = std(dumT);
  end
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
 
  end % end loop_ens
  
  EnRStat(1,loop_lam) = mean(EnR);
  EnRStat(2,loop_lam) = std(EnR);
  EnTStat(1,loop_lam) = mean(EnT);
  EnTStat(2,loop_lam) = std(EnT);
  
 end % end loop_lam
  
 if COMM
 cprintf([0.3,0.3,0.3],'-----------    END: Combine_Rows   ------------\n')
 cprintf([0.3,0.3,0.3],'-----------------------------------------------\n')
 end
end % end if



return