% fn_AttnModels
%
% INPUTS:
%
% Hm = scalar/vector wave heights [mm]
% Tp = scalar/vector periods      [s]
%
% FLAGS:
%
% ATTN = plot (dimensional) attenuation coefficient?
% EN   = plot transmitted energy?
% DIRS = plot transmitted energy as function of angle?
% EVALS = plot eigenvalues of periodic problem?
% PRBS = [3d real geom, 3d random phases, 2d, Boltzmann steady]
% TAMP = Transmitted amplitudes
% LONG = long floe limit (for 2d code)
% COMM = comments on (1) or off (0)

function out=Main_AttnModels(Tp,conc,PRBS,Nd,COMM,DO_PLOT)

%%%%%%%%%%%%%%%%%%%%%%
%% %%%% PRELIMS %%%%%%
%%%%%%%%%%%%%%%%%%%%%%

%% GENERAL:

if ~exist('PRBS','var');    PRBS='Boltzmann steady'; end
if ~exist('COMM','var');    COMM=1;         end
if ~exist('DO_PLOT','var'); DO_PLOT=0;      end

if ~exist('Tp','var'); Tp = .65; end
if ~exist('conc','var'); conc=79; end

if ~exist('rigid','var'); rigid=4; end
 
if ~exist('Nd','var'); Nd = 1e2; end

if ~exist('th_res','var'); th_res=200; end

if ~exist('TEST','var'); TEST='Oceanide'; end

if ~exist('SURGE','var'); SURGE=1; end

if ~exist('extra_pts','var'); extra_pts=[]; end
if ~exist('terms_grn','var'); terms_grn=100; end

if strcmp(TEST,'Oceanide')
    
 if ~exist('Param','var'); Param = ParamDef_Oceanide(rigid); 
    Param = ModParam_def(Param,1,Nd,extra_pts,terms_grn,th_res); end
 
else
    
 cprintf('red',['Not set up yet' '\n'])
 
end

conc=conc/100;
out_str = ' ''Dummy'' ';
out_val = ' 0 ';

%% 3D PROBLEM:

if conc==.79
 file_mkr0 = 'ela_plt_Rigid16epspt005_T';
else
 file_mkr0 = 'ela_plt_Rigid8epspt1_T';
end

%% 2D PROBLEM:

%% BOLTZMANN:

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%% 2d attenuation (long floe limit & no long floe limit) %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Long floe limit

if strfind(PRBS,'2d lfl EMM')
 
 attn_2d = zeros(1,length(Tp));
 
 for loop_p=1:length(Tp)
  out = fn_ElasticRaft2d('freq',1/Tp(loop_p),Param,'transmitted energy',...
   SURGE,1,COMM);
  attn_2d(loop_p) = out(1).value; clear out
 end
 
 T_2dl = exp(conc*log(attn_2d)*Param.MIZ_length/Param.floe_diam/2); 
 clear attn_2d
 
 out_str = [out_str '; ''2d long'' '];
 out_val = [out_val '; T_2dl '];
 
end % end if 2d long

%%% No long floe limit

if strfind(PRBS,'2d EMM') 
 
 if ~exist('ens1','var'); ens1 = 1; end
 
 attn_2d = zeros(1,length(Tp));
 
 for loop_p=1:length(Tp)
   out = fn_ElasticRaft2d('freq',1/Tp(loop_p),Param,'transmitted energy',...
    SURGE,0,COMM,ens1);
   attn_2d(loop_p) = out(1).value; clear out
 end
  
 T_2dx = exp(conc*log(attn_2d)*Param.MIZ_length/Param.floe_diam/2);  
 clear attn_2d
 
 out_str = [out_str '; ''2d no long'' '];
 out_val = [out_val '; T_2dx '];
 
end % end if 2d no long

%%% Integral equation method (from file): 

if strfind(PRBS,'2d BIE')
 
 if strcmp(getenv('LOGNAME'),'a1612881')
  basedir='/Volumes/scratch/Xtra_Data/CITEPH-64-2012/Model/SingleFloe/';
 end
 disp('2d model:')
 dum_dirs = fn_FolderNames(basedir);
 for loop_d=1:length(dum_dirs)
  disp([int2str(loop_d) '. ' dum_dirs{loop_d}])
 end
 ind = str2num(input('Choose folder number... ', 's'));
 basedir=[basedir dum_dirs{ind} '/'];
 dum_dirs = fn_FileNames(basedir);
 for loop_d=1:length(dum_dirs)
  disp([int2str(loop_d) '. ' dum_dirs{loop_d}])
 end
 cprintf(0.4*[1,1,1],'xxxxxxxxxxx-Xint-Xmodes-Nfac\n')
 ind = str2num(input('Choose file number... ', 's'));
 basedir=[basedir dum_dirs{ind}];
 load(basedir,'out','for_vec','v_info')
 cprintf(0.4*[1,1,1],[v_info{1}(1:end-1) '\n'])
 clear ind basedir dum_dirs
 count=1;
 while and(count>0,count<=length(out))
  if strfind(out(count).name,'reflected energy')
   attn_2d=-log(1-out(count).value);
   count=0;
  else
   count=count+1;
  end % END IF
 end % END WHILE
 
 T_2di = exp(-conc*attn_2d*Param.MIZ_length/Param.floe_diam/2); clear attn_2d
 
 out_str = [out_str '; ''2d BIE'' ; ''periods'' '];
 out_val = [out_val '; T_2di ; for_vec'];
 
 clear v_info count out
  
%  for loop_p=1:length(Tp)
%   out = fn_Solve_2dTEP_Cos('freq',1/Tp(loop_p),Param,0,COMM,0);
%   for loop_out=1:length(out)
%    if strcmp(out(loop_out).name,'Amp_out')
%     attn_2d(loop_p) = -log(1-abs(out(loop_out).value)^2);
%    end
%   end
%  end
%  
%  attn_2d = conc*attn_2d/Param.floe_diam/2;
%  
%  T_2di = exp(-attn_2d*Param.MIZ_length); clear attn_2d
%  
%  out_str = [out_str '; ''2d IE'' '];
%  out_val = [out_val '; T_2di '];
 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%   3d conceptual     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if strfind(PRBS,'3d simple')
 
 if ~exist('th_vec','var')
  th_vecA=linspace(0,pi/2,101);
  th_vec=unique([th_vecA,th_vecA+pi/2]); th_vec=unique([th_vec,-th_vec]);
  clear th_vecA
  indsT=find(abs(th_vec)<pi/2); indsR=find(or(th_vec<-pi/2,th_vec>pi/2));
%   th_vecR=th_vecA+pi/2; th_vecR=[th_vecR,nan,-fliplr(th_vecR)];
%   th_vecT=unique([th_vecA,-th_vecA]); clear th_vecA
 end
 
 for loop_p=1:length(Tp)
  out = fn_ElasticDisk('freq',1/Tp(loop_p), Param, 'DTM', ...
   th_vec, 0, SURGE, COMM, 0);
  count=1;
  while and(count>0,count<=length(out))
  if strfind(out(count).name,'DTM')
   DTM=out(count).value;
   count=0;
  else
   count=count+1;
  end % END IF
 end % END WHILE
 
 M = (length(DTM)-1)/2; M_vec=-M:M; 
 Ip = exp(1i*M_vec*pi/2);
 Inc=0*th_vec;
 for loop_th=1:length(th_vec)
  Inc(loop_th) = sum(Ip.*exp(1i*M_vec*(th_vec(loop_th)-pi/2)))/2/pi;
 end
 Hp = exp(-1i*M_vec*pi/2);
 Rsc=0*th_vec(indsR); Tsc=0*th_vec(indsT);
 ct_th=1;
 for loop_th=indsR
  Rsc(ct_th) = sum(Ip.*Hp.*DTM.*exp(1i*M_vec*th_vec(loop_th)))/pi; 
  ct_th=ct_th+1; 
 end
 ct_th=1;
 for loop_th=indsT
  Tsc(ct_th) = sum(Ip.*Hp.*DTM.*exp(1i*M_vec*th_vec(loop_th)))/pi; 
  ct_th=ct_th+1;
 end
 I = fn_TrapRule(th_vec,abs(cos(th_vec)).*abs(Inc).^2);
 [thR,inds] = unique(mod(th_vec(indsR),2*pi)); 
 R = fn_TrapRule(thR,abs(cos(thR)).*abs(Rsc(inds)+Inc(indsR(inds))).^2);
 T = fn_TrapRule(th_vec(indsT),abs(cos(th_vec(indsT))).*abs(Tsc+Inc(indsT)).^2);
 clear DTM count M M_vec Ip
 if COMM
  cprintf('blue',['R/I = ' num2str(R/I) '\n'])
  cprintf('blue',['T/I = ' num2str(T/I) '\n'])
  if abs(1-(R+T)/I)>1e-3
   cprintf('green',['energy error: Tp=' num2str(Tp(loop_p)) ...
    ' (' num2str(abs(1-(R+T)/I)) ')\n'])
  end
 end % end if COMM
 end % END FOR LOOP_P
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%      Boltzmann steady      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if strfind(PRBS,'Boltzmann steady')
 
 if ~exist('th_res','var'); th_res=1*50; end
 
 %I_Boltz = zeros(length(Tp),th_res);
 
 for loop_p=1:length(Tp) 
  out = ...
   fn_Boltzmann_Steady('freq',1/Tp(loop_p),conc,Param,...
   'int trans energy',COMM,0);
  T_Blt(loop_p) = out(1).value; clear out
 end
 
%  [~,ind] = min(abs(th_vec));
%  
%  T_Blt = sqrt(I_Boltz(:,ind));

 T_Blt = sqrt(T_Blt);
 
 out_str = [out_str '; ''Boltzmann steady'' '];
 out_val = [out_val '; T_Blt '];
 
end % end if Boltzmann steady

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%   3d Wavetank (real geometry)   %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if strfind(PRBS,'3d WT')
 
 if ~exist('Rows','var'); Rows = 5; end
 if ~exist('ens0','var'); ens0 = 100; end
 
 EnTStat = zeros(6,length(Tp)); attn_3d = zeros(6,length(Tp));
 
 for loop_p=1:length(Tp)
  
  file_marker = [file_mkr0 int2str(floor(Tp(loop_p))) 'pt' ...
   int2str(10*(Tp(loop_p)-floor(Tp(loop_p))))];
  
  [TraStat{loop_p},EnTStat(:,loop_p),attn_3d(:,loop_p)]=...
   fn_CombineRows(file_marker,Rows,ens0,0,COMM);
  
  M(loop_p) = length(TraStat{loop_p}(1,:));
  
 end
 
 out_str = [out_str '; ''3d WT'' '];
 out_val = [out_val '; T_3d '];
 
end % end 3d prb

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%% 3d Wavetank (random phases between rows)  %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if strfind(PRBS,'3d rnd WT')
 
 EnTStat = zeros(3,length(Tp)); attn_3d = zeros(4,length(Tp));
 
 for loop_p=1:length(Tp)
  
  file_marker = [file_mkr0 int2str(floor(Tp(loop_p))) 'pt' ...
   int2str(10*(Tp(loop_p)-floor(Tp(loop_p))))];
  
  [TraStat{loop_p},EnTStat(:,loop_p),attn_3d(:,loop_p)]=...
   fn_CombineRows(file_marker,Rows,ens0,10*pi/Tp(end),COMM);
  
 end
 
 %%% PLOTS
 
 %%%% Attenuation coeff
 
 if ATTN
  figure(ATTN)
  plot(Tp,attn_3d(1,:),[col_3db mkr],'markersize',12)
  plot(Tp,attn_3d(2,:),[col_3db '*'],'markersize',10)
  plot(Tp,attn_3d(3,:),[col_3db '.'],'markersize',10)
  plot(Tp,attn_3d(4,:),[col_3db 'x'],'markersize',10)
 end
 
 %%%% Transmitted energy
 
 if EN
  figure(EN)
  plot(Tp,EnTStat(1,:),[col_3db mkr],'markersize',12)
 end
 
 %%%% Directional spectrum
 
 if DIRS
  for loop_p=1:length(Tp)
   figure(DIRS+loop_p)
   plot(TraStat{loop_p}(3,:),TraStat{loop_p}(1,:),[col_3db mkr],'markersize',8)
   plot(TraStat{loop_p}(3,:),TraStat{loop_p}(3,:),[col_3db '*'],'markersize',8)
  end
 end
 
 out_str = [out_str '; ''3d rand WT'' '];
 out_val = [out_val '; T_3dr '];
 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%% Rows (random phases between rows)  %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if strfind(PRBS,'Rows')
 
 rws = 1;
 
 for loop_p=1:length(Tp) 
  out = fn_Rows(Param,'freq',1/Tp(loop_p),conc,rws,1,'transmitted energy',COMM);
  attn_rows(loop_p) = out(1).value; clear out   
 end
 
 T_row = exp(conc*log(attn_rows)*Param.MIZ_length/Param.floe_diam/2);  
 clear attn_2d
  
 out_str = [out_str '; ''Rows'' '];
 out_val = [out_val '; T_row '];
 
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%            Output          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

eval(['out=struct( ''name'', {' out_str ...
 '}, ''value'', {' out_val '});'])

out(1)=[];

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SUB-FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function I = fn_TrapRule(xx,yy)

xtil=diff(xx);
ytil=0.5*(yy(1:end-1)+yy(2:end));

I = sum(xtil.*ytil);

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% if ~exist('col','var');     mkr='o'; end
% if ~exist('col_2d','var');  col_2d='b'; end
% if ~exist('col_3d','var');  col_3d='r'; end
% if ~exist('col_3db','var'); col_3db='g'; end
% if ~exist('col_Bol','var'); col_Bol='m'; end
% 
% if DO_PLOT
%  if ~exist('ATTN','var');  ATTN=0*101;  end
%  if ~exist('EN','var');    EN=0*102;    end
%  if ~exist('DIRS','var');  DIRS=0*300;  end
%  if ~exist('EVALS','var'); EVALS=0*401; end
%  if ~exist('TAMP','var');  TAMP=1*501; end
% else
%  ATTN=0; EN=0; DIRS=0; EVALS=0; TAMP=0;
% end
% 
% if ATTN;  figure(ATTN);  hold on; end
% if TAMP;  figure(TAMP);  hold on; end
% if EN;    figure(EN);    hold on; end
% if EVALS; figure(EVALS); hold on; end
% if DIRS;
%  for loop_p=1:length(Tp)
%   figure(DIRS+loop_p); hold on;
%   title(['period = ' num2str(Tp(loop_p))])
%   set(gca,'yscale','log');
%   set(gca,'ylim',[0,1])
%  end
% end
% 
% %%% PLOTS
% 
% %%%% Attenuation coeff
% 
% if ATTN
%  figure(ATTN)
%  plot(Tp,attn_2d,'bo','markersize',12)
% end
% 
% %%%% Transmitted amplitude
% 
% if TAMP
%  figure(TAMP)
%  plot(Tp,(Hm/2).*exp(-attn_2d*5),'bo','markersize',12)
% end
% 
% %%% PLOTS
% 
% %%%% Attenuation coeff
% 
% if ATTN
%  figure(ATTN)
%  plot(Tp,attn_2d,'rx','markersize',12)
% end
% 
% %%%% Transmitted amplitude
% 
% if TAMP
%  figure(TAMP)
%  plot(Tp,(Hm/2).*exp(-attn_2d*5),'rx','markersize',12)
% end
% 
% if DIRS
%  for loop_p=1:length(Tp)
%   figure(DIRS+loop_p)
%   plot(cos(th_vec),I_Boltz(loop_p,:),[col_Bol '-.'])
%  end
%  %  for loop_p=2:length(pers)
%  %   plot(cos(th_vec),(loop_p-1)+0*I_Boltz(loop_p,:),'k:')
%  %  end
%  %  set(gca,'box','on')
% end
% 
% %%% PLOTS
% 
% %%%% Attenuation coeff
% 
% if ATTN
%  figure(ATTN)
%  plot(Tp,attn_3d(1,:),[col_3d mkr],'markersize',12)
%  plot(Tp,attn_3d(2,:),[col_3d '*'],'markersize',10)
%  plot(Tp,attn_3d(3,:),[col_3d '.'],'markersize',10)
%  plot(Tp,attn_3d(4,:),[col_3d 'x'],'markersize',10)
%  title(['full perturbed (' mkr '); periodic (*); periodic eval (.); zeroth order (x)'])
% end
% 
% %%%% Transmitted energy
% 
% if EN
%  figure(EN)
%  plot(Tp,EnTStat(1,:),[col_3d mkr],'markersize',12)
% end
% 
% %%%% Directional spectrum
% 
% if DIRS
%  for loop_p=1:length(Tp)
%   figure(DIRS+loop_p)
%   plot(TraStat{loop_p}(7,:),TraStat{loop_p}(1,:),['bo'],'markersize',8)
%   plot(TraStat{loop_p}(7,:),TraStat{loop_p}(5,:),['bx'],'markersize',8)
%   plot(TraStat{loop_p}(7,:),TraStat{loop_p}(3,:),['rs'],'markersize',8)
%   plot(TraStat{loop_p}(7,:),TraStat{loop_p}(4,:),['r+'],'markersize',8)
%  end
% end
% 
% %%%% eigenvalues of periodic problem
% 
% if EVALS
%  figure(EVALS)
%  for loop_p=1:length(Tp)
%   dum(loop_p) = max(abs(imag(TraStat{loop_p}(9,:))));
%  end
%  plot(Tp,dum,[col_3d mkr],'markersize',12); clear dum
% end