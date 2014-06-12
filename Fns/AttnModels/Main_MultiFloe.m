% function Main_MultiFloe(Np,fortyp,lam_vec,scatyp,row_seps,file_marker,COMM,POOL)
%
% Calculates and save reflection/transmission by single rows of scatterers
% enclosed by walls.
% Location of scatterers is randomly chosen.
% For later use in <enter file here later>
%
% INPUTS:
%
% TEST = sting identifying the test to be performed
%        options: `Oc79', `Oc39', 'Quick'
% Np = number of floes
% Ndtm = no. vert modes used in calc DTMs
% fortyp = 'freq' or 'wlength' or 'waveno'
% lam_vec = a vector of forcing fortyp
% scatyp = 'none' (no scatterers) or 'soft' (sound soft cylinders) or 'cyl'
%           (sound hard cylinders) or `ela_plt' (elastic plates)
% ens = ensemble size
% epsilon = perturbation from periodic position (dimensional)
% file_marker = string for file identifier; use 0 for no write
% COMM = flag for comments on (1) or off (0)
% POOL = flag for matlabpool on (number of cpus) or off (0)

function Main_MultiFloe(TEST,Ndtm,fortyp,lam_vec,ens,file_marker,add_per,COMM,POOL)

if ~exist('POOL','var'); POOL=0; end
if ~exist('TEST','var'); TEST='Oc39'; end
if ~exist('COMM','var'); COMM=1; end

if ~exist('extra_pts','var'); extra_pts=[]; end
if ~exist('terms_grn','var'); terms_grn=100; end
if ~exist('fortyp','var'); fortyp='freq'; end
if ~exist('lam_vec','var'); lam_vec=1/2; end 
if ~exist('Ndtm','var'); Ndtm=5; end

if ~exist('add_per','var'); add_per=0; end

if length(lam_vec)>1
 disp('DONT SOLVE MORE THAN 1 PROBLEM AT A TIME !!!!')
 return
end

if POOL
 matlabpool('open',POOL)
end

if ~exist('file_marker','var'); file_marker='test'; end

if strcmp(TEST,'Oc79')
 
 if ~exist('RIGID','var'); rigid=10; end
 if ~exist('scatyp','var'); scatyp='ela_plt'; end
 if ~exist('Np','var'); Np=5*16; end
 if ~exist('ens','var'); ens=1; end
 if ~exist('scat_ident','var'); scat_ident=1; end
 if strcmp(fortyp,'freq')
  if ~exist('epsilon','var'); epsilon=fn_surge(1/lam_vec); end
 end
 str_GeomDisks = 'fn_GeomDisks_Oc79';
 Param = ParamDef_Oceanide(rigid,Np);
 Param = ModParam_def(Param,1,Ndtm,extra_pts,terms_grn);
 TankDim = [40 16 Param.bed];
 
elseif strcmp(TEST,'Oc39')
 
 if ~exist('RIGID','var'); rigid=10; end
 if ~exist('scatyp','var'); scatyp='ela_plt'; end
 if ~exist('Np','var'); Np=5*8; end
 if ~exist('ens','var'); ens=1; end
 if ~exist('scat_ident','var'); scat_ident=1; end
 if strcmp(fortyp,'freq')
  if ~exist('epsilon','var'); epsilon=fn_surge(1/lam_vec); end
 end
 str_GeomDisks = 'fn_GeomDisks_Oc39';
 Param = ParamDef_Oceanide(rigid,Np);
 Param = ModParam_def(Param,1,Ndtm,extra_pts,terms_grn);
 TankDim = [40 16 Param.bed]; 
 
elseif strcmp(TEST,'Quick')
 
 if ~exist('RIGID','var'); RIGID=1; end
 if ~exist('scatyp','var'); scatyp='ela_plt'; end
 if ~exist('Np','var'); Np=2; end
 if ~exist('ens','var'); ens=1; end
 if ~exist('epsilon','var'); epsilon=0.5; end
 if ~exist('scat_ident','var'); scat_ident=1; end
 str_GeomDisks = 'fn_GeomDisks_Quicky';
 TankDim = [40 16 3]; 
 
else
    
 disp('code something else!!!')
 
end

disp('%--------------------------------------------%')
disp('%---------- START: Main_MultiFloe -----------%')

disp(['Test: ' TEST])
if strcmp(fortyp,'freq')
    disp(['freq = ' num2str(lam_vec) ' Hz'])
end
disp([int2str(ens) ' in ensemble'])
%clockout=clock;
% display(['Started: ' int2str(clockout(3)),'/', int2str(clockout(2)),'/',int2str(clockout(1)),' ' ,...
%   int2str(clockout(4)),':', int2str(clockout(5)),':',int2str(round(clockout(6)))] );
%clear clockout

disp('Started:') ; 
fn_dispclock;

tic

%% Creation of local paths
%CreatePaths_Mac

%% Wavelengths
%%% NB. There are m `tank modes' when k\in \pi(m,m+1)/width

%% Parametrisation
% Note that the origin of the Cartesian system is located at the bottom 
% left corner of the wave tank and coincides with the free-surface at rest.

%%% Geometry of the wave tank: [length width height] %%%


%%% Unperturbed problem

% eval(['GeomDisks=' str_GeomDisks '(TankDim,0,COMM,100);']);
%  
% %%% Physical properties
% Param = ParamDef3d(GeomDisks,RIGID);
% 
% %%% Generation of horizontal mesh points for the free-surface region
% Mesh = Mesh_FS_def(GeomDisks, TankDim);
% 
% %%% Modal parameters and accuracy
% Param = ModParam_def(Param,1,Ndtm,extra_pts,terms_grn);
% 
% [Rm_sv0,Tm_sv0,Rp_sv0,Tp_sv0]...
%     = fn_MultiMode_MultiFloe_loop(lam_vec, Param, TankDim, ...
%     GeomDisks, Mesh, scatyp, fortyp, COMM);
% 
% GeomDisks_sv0 = GeomDisks;

for loop_ens=1:ens

 if add_per
  if loop_ens~=1
   eps0=epsilon;
  else
   eps0=0;
  end
 else
  eps0=epsilon;
 end

if strcmp(TEST,'Oc79')
 GeomDisks=fn_GeomDisks_Oc79(TankDim,eps0,COMM,100);
elseif strcmp(TEST,'Oc39')
 GeomDisks=fn_GeomDisks_Oc39(TankDim,eps0,COMM,100);
elseif strcmp(TEST,'Quick')
 GeomDisks=fn_GeomDisks_Quicky(TankDim,eps0,COMM,100);
 Param = ParamDef3d(GeomDisks,RIGID);
 Param = ModParam_def(Param,1,Ndtm,extra_pts,terms_grn);
end

%%% Physical properties

%%% Generation of horizontal mesh points for the free-surface region
%Mesh = Mesh_FS_def(GeomDisks, TankDim);

%%% Modal parameters and accuracy


[Rm_sv{loop_ens},Tm_sv{loop_ens},Rp_sv{loop_ens},Tp_sv{loop_ens},...
      v_sv{loop_ens},k0_sv{loop_ens},x_lim_sv{loop_ens},reson_sv{loop_ens}]...
    = fn_MultiMode_MultiFloe_loop(lam_vec, Param, TankDim, ...
    GeomDisks, scatyp, fortyp, scat_ident, COMM);

GeomDisks_sv{loop_ens} = GeomDisks;

end % end loop_ens

if ens>1
 k0_sv = k0_sv{1}; v_sv = v_sv{1}; reson_sv = reson_sv{1}; 
 for loop=1:ens
  x_lim_sv{loop} = x_lim_sv{loop}{1};
 end
end

tm=toc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%         SAVE           %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

if file_marker~=0
    
 %%% CHANGE TO LOCAL DIRECTORY' %%%   
 %file_pre = '../../../../../Documents/MatLab/Data/Wavetank/Full/';
 file_pre = '../../../Data/CITEPH-64-2012/MultiFloes/';
 if strcmp(TEST,'Oc79')
  file_pre=[file_pre 'Conc79/'];
 elseif strcmp(TEST,'Oc39')
  file_pre=[file_pre 'Conc39/'];
 elseif strcmp(TEST,'Quick')  
  file_pre=[file_pre 'Quick/'];
 end
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
 v_info = {'Version 1.1:'; 
     ['31.10.2013: variation all in x direction'];
     ['17.07.2013: modified from Main_Rows']};
     
 file_name = ['Main_MultiFloe_',scatyp,'_',file_marker,'.mat'];
   
 eval(['GeomDisks=' str_GeomDisks '(TankDim,0,COMM,100)']);
 Param = ParamDef3d(GeomDisks);
 Mesh = Mesh_FS_def(GeomDisks, TankDim);
 Param = ModParam_def(Param,extra_pts,terms_grn);
 
 hrs = floor(tm/3600);mins = floor((tm - hrs*3600)/60);secs = tm - hrs*3600 - mins*60;
 clockout = clock; duration = [hrs,mins,secs];
    
 for loop_ens=1:ens
  for loopi=1:4
   dum_str{loopi} = num2str(GeomDisks_sv{loop_ens}(1,loopi));    
   for loopii=2:Np   
    dum_str{loopi} = [ dum_str{loopi}, ', ' ...
        num2str(GeomDisks_sv{loop_ens}(loopii,loopi)) ];
   end
  end
  clear loopi loopii
  dum_str{5} = num2str(Param.E(1));    
  for loopii=2:Np   
   dum_str{5} = [ dum_str{5}, ', ' num2str(Param.E(loopii)) ];
  end
  clear loopii
  disks{loop_ens} = ['x=' dum_str{1} ...
     '; y=' dum_str{2} ...
     '; r=' dum_str{3} ...
     '; h=' dum_str{4} ...
     '; E=' dum_str{5} ];
 end % end loop_ens
 clear dum_str
  
 description = {['Problem: ' scatyp '; ensemble size = ' int2str(ens)], ...
     [fortyp ' =' num2str(lam_vec(1))], ...
     ['eps = ', num2str(epsilon)] ... 
     ['Tank: ' num2str(TankDim(1)) 'x' num2str(TankDim(2)) ...
     'x' num2str(TankDim(3))], ...
     [int2str(Np) ' Disk(s)' ], ...
     ['Vert modes (int & DTM)=' int2str(Param.Nint) ...
     ' & ' int2str(Param.Ndtm) ...
     ', extra points=' int2str(Param.extra_pts) ...
     ', angular res=' int2str(Param.res_green) ...
     ', Green fn terms=' int2str(Param.terms_green) ...
     ', Greens fn cutoff=' int2str(Param.cutoff_green) ...
     ', resonance tol=' int2str(Param.tolres)]};
 
%  [fortyp ' :' num2str(lam_vec(1)) ' -> ' num2str(lam_vec(end)) ' (' ...
%       int2str(length(lam_vec)) ')'] ...

 w = TankDim(2);
 
 save([file_pre,file_name], 'v_info', 'lam_vec', 'fortyp', ...
    'Rm_sv', 'Tm_sv', 'Rp_sv', 'Tp_sv', 'v_sv', 'w', ...
    'description', 'x_lim_sv', 'epsilon', 'disks', 'k0_sv', ...
    'reson_sv', 'GeomDisks_sv', 'clockout', 'duration')
%,... 'Rm_sv0','Tm_sv0','Rp_sv0','Tp_sv0','GeomDisks_sv0')

else
 beep; beep; beep
end

if POOL; matlabpool close; end

if COMM
disp('%----------- END: Main_MultiFloe ------------%')
disp('%--------------------------------------------%')
end

return

%%%%%%%%%%%%%%%%%%%%
%%% SUBFUNCTIONS %%%
%%%%%%%%%%%%%%%%%%%%

function [Rm,Tm,Rp,Tp,v_vec,k0,x_lim,reson_mkr]...
    = fn_MultiMode_MultiFloe_loop(lam_vec, Param, TankDim, ...
    GeomDisks, scatyp, fortyp, scat_ident, COMM)

for loop_lam=1:length(lam_vec)
         
 %%% Forcing
 Forcing = Force_def(Param.g(1), TankDim(3), fortyp, lam_vec(loop_lam));

 %% Solution - Disks in a channel
 % Scattering matrix is S = [Rm Tp;Tm Rp]
 % Wavenumbers matrices v_vec and u_vec are scaled by vertical wavenumbers.
 [Rm{loop_lam},Tm{loop_lam},Rp{loop_lam},Tp{loop_lam},...
     v_vec{loop_lam},~,k0{loop_lam},~,x_lim{loop_lam},reson_mkr{loop_lam}] = ...
    fn_MultiMode_MultiFloe([Param.rho_0 Param.rho Param.nu Param.E ...
    Param.g], Param.Nint, Param.Ndtm, Param.Mev, TankDim, Forcing.kappa, ...
    Param.beta.', GeomDisks(:,4).', Param.draft.', GeomDisks(:,3).', ...
    [GeomDisks(:,1).'; GeomDisks(:,2).'], ...
    [Param.res_green, Param.terms_green, Param.cutoff_green, Param.tolres], ...
    Param.extra_pts, scatyp, scat_ident, COMM);
   %Mesh.r_vec, Mesh.th_vec, Mesh.x_vec, Mesh.y_vec, Mesh.FS_mesh, ...
    
 end % end loop_lam
 
return

%
 
function GeomDisks=fn_GeomDisks_Oc79(TankDim,epsilon,COMM,mxc)

count=1; flg=1;
Np = 16;

while and(flg~=0,count<=mxc)
 %%% Location and geometry of disks: [x_c y_c Rad thick]
 for loop_rows=1:5
  for loop=1:Np
   epsx = 2*(rand-0.5)*epsilon; epsy = 2*(rand-0.5)*epsilon;
   GeomDisks(Np*(loop_rows-1)+loop,:) = ...
      [17.5+(loop_rows-1)+epsx epsy+(TankDim(2)/2/Np)+(TankDim(2)/Np)*(loop-1) 0.495 33e-3];
  end
 end

 %%% Check disk/disk and disk/wall overlaps
 flg=overlap(GeomDisks, TankDim, COMM);
end

if count>mxc
 disp('Cannot get plates configured without overlap!!!')
 disp('Exiting')
 exit
end

return

%

function GeomDisks=fn_GeomDisks_Oc39(TankDim,epsilon,COMM,mxc)

count=1; flg=1;
Np = 8;

while and(flg~=0,count<=mxc)
 %%% Location and geometry of disks: [x_c y_c Rad thick]
 for loop_rows=1:2:5
  for loop=1:Np
   epsx = 2*(rand-0.5)*epsilon; epsy = 0; %2*(rand-0.5)*epsilon;
   GeomDisks(Np*(loop_rows-1)+loop,:) = ...
      [17.5+(loop_rows-1)+epsx epsy+(TankDim(2)/2/2/Np)+(TankDim(2)/2/Np)*(2*loop-1) 0.495 33e-3];
  end
 end
 
 for loop_rows=2:2:4
  for loop=1:Np
   epsx = 2*(rand-0.5)*epsilon; epsy = 0; %2*(rand-0.5)*epsilon;
   GeomDisks(Np*(loop_rows-1)+loop,:) = ...
      [17.5+(loop_rows-1)+epsx epsy+(TankDim(2)/2/2/Np)+(TankDim(2)/2/Np)*(2*loop-2) 0.495 33e-3];
  end
 end

 %%% Check disk/disk and disk/wall overlaps
 flg=overlap(GeomDisks, TankDim, COMM);
end

if count>mxc
 disp('Cannot get plates configured without overlap!!!')
 disp('Exiting')
 exit
end

return

% 

function GeomDisks=fn_GeomDisks_Quicky(TankDim,epsilon,COMM,mxc)

count=1; flg=1;
Np = 2;

while and(flg~=0,count<=mxc)
 %%% Location and geometry of disks: [x_c y_c Rad thick]
 for loop=1:Np
   epsx = 2*(rand-0.5)*epsilon; epsy = 2*(rand-0.5)*epsilon;
   GeomDisks(loop,:) = ...
      [17.5+epsx epsy+(TankDim(2)/2/Np)+(TankDim(2)/Np)*(loop-1) 0.495 33e-3];
 end

 %%% Check disk/disk and disk/wall overlaps
 flg=overlap(GeomDisks, TankDim, COMM);
end

if count>mxc
 disp('Cannot get plates configured without overlap!!!')
 disp('Exiting')
 exit
end

return

% 

function sg=fn_surge(tau)

pers = [0.65,0.95,1.25,1.55,1.85];

%expt_sg = [0.12,0.50,0.86,0.98,1.03];

expt_sg = [0.08,1.17,3.06,3.81,3.98]./100;

if tau<0.65
 sg=0;
elseif tau>1.85
 sg=3.98/100;
else
 sg=interp1(pers,expt_sg,tau);
end

return