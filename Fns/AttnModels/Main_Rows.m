% function Main_Rows(Np,fortyp,lam_vec,scatyp,row_seps,file_marker,COMM,POOL)
%
% Calculates and save reflection/transmission by single rows of scatterers
% enclosed by walls.
% Location of scatterers is randomly chosen.
% For later use in <enter file here later>
%
% INPUTS:
%
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

function Main_Rows(Np,Ndtm,fortyp,lam_vec,scatyp,conc,file_marker,COMM,POOL)

if ~exist('POOL','var'); POOL=0; end
if ~exist('COMM','var'); COMM=1; end
if ~exist('RIGID','var'); RIGID=4; end

if POOL
 matlabpool('open',POOL)
end

if ~exist('Np','var'); Np=16; end
if ~exist('extra_pts','var'); extra_pts=[]; end
if ~exist('terms_grn','var'); terms_grn=100; end
if ~exist('file_marker','var'); file_marker='test'; end
if ~exist('fortyp','var'); fortyp='freq'; end
if ~exist('ens','var'); ens=0; end
if ~exist('epsilon','var'); epsilon=0; end
if ~exist('lam_vec','var'); lam_vec=1/2; end %4.0317;1.0472;%2*pi-0.1:0.1:2*pi+0.1;0.95*pi/16:0.05*pi/16:1.05*pi/16;
if ~exist('scatyp','var'); scatyp='ela_plt'; end
if ~exist('Ndtm','var'); Ndtm=5; end

if length(lam_vec)>1
 disp('DONT SOLVE MORE THAN 1 PROBLEM AT A TIME !!!!')
 exit
end

if conc==79
 w=1;
end

if ~exist('w','var'); w=16; end

disp('%---------------------------------------%')
disp('%---------- START: Main_Rows -----------%')

disp([int2str(Np) ' Scatterer(s)'])
if strcmp(fortyp,'freq')
    disp(['freq = ' num2str(lam_vec) ' Hz'])
end
disp([int2str(ens) ' in ensemble'])

tic

%% Creation of local paths
%CreatePaths_Mac

%% Wavelengths
%%% NB. There are m `tank modes' when k\in \pi(m,m+1)/width

%% Parametrisation
% Note that the origin of the Cartesian system is located at the bottom 
% left corner of the wave tank and coincides with the free-surface at rest.

%%% Geometry of the wave tank: [length width height] %%%
cell_len=w/Np; % side length of `cell' containing a floe = width/no. floes
TankDim = [cell_len w 3.1]; %[15 10 0.5]; % [40 16 2];

%%% Unperturbed problem

GeomDisks=fn_GeomDisks(Np,TankDim,0,COMM,100);
 
%%% Physical properties
%Param = ParamDef3d(GeomDisks,RIGID);
%%% Modal parameters and accuracy
%Param = ModParam_def(Param,1,Ndtm,extra_pts,terms_grn);

if ~exist('Param','var'); Param = ParamDef_Oceanide(RIGID); 
    Param = ModParam_def(Param,1,Ndtm,extra_pts,terms_grn); end

%%% Generation of horizontal mesh points for the free-surface region
Mesh = Mesh_FS_def(GeomDisks, TankDim);

[Rm_sv0,Tm_sv0,Rp_sv0,Tp_sv0]...
    = fn_MultiMode_MultiFloe_loop(lam_vec, Param, TankDim, ...
    GeomDisks, Mesh, scatyp, fortyp, COMM);

GeomDisks_sv0 = GeomDisks;

for loop_ens=1:ens

GeomDisks=fn_GeomDisks(Np,TankDim,epsilon,COMM,100);
 
%%% Physical properties
Param = ParamDef3d(GeomDisks,RIGID);

%%% Generation of horizontal mesh points for the free-surface region
Mesh = Mesh_FS_def(GeomDisks, TankDim);

%%% Modal parameters and accuracy
Param = ModParam_def(Param,1,Ndtm,extra_pts,terms_grn);

[Rm_sv{loop_ens},Tm_sv{loop_ens},Rp_sv{loop_ens},Tp_sv{loop_ens},...
      v_sv{loop_ens},k0_sv{loop_ens},x_lim_sv{loop_ens},reson_sv{loop_ens}]...
    = fn_MultiMode_MultiFloe_loop(lam_vec, Param, TankDim, ...
    GeomDisks, Mesh, scatyp, fortyp, COMM);

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
 if conc==79
  file_pre = '/home/users/lbennetts/Data/CITEPH-64-2012/SingleFloe/Conc79';
 end
 %file_pre = '../../Data/Wavetank/Rows/';
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
 v_info = {'Version 1.6:'; 
     ['20.01.2014: quasi-periodic problem'];
     ['17.07.2013: unperturbed problem added'];
     ['16.07.2013: thickness now 33e-3 not 33e-2'];
     ['09.07.2013: radius now .495 not .5'];
     ['19.06.2013: separated vert dims for interactions & DTMs'];
     ['11.06.2013: restricted to one frequency at a time'];
     ['22.05.2013: modified from Main_Channel_freq']};
     
 file_name = ['Main_Rows_',scatyp,'_',file_marker,'.mat'];
   
 GeomDisks=fn_GeomDisks(Np,TankDim,0,COMM,100);
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
 
 save([file_pre,file_name], 'v_info', 'lam_vec', 'fortyp', ...
    'Rm_sv', 'Tm_sv', 'Rp_sv', 'Tp_sv', 'v_sv', 'w', ...
    'description', 'x_lim_sv', 'epsilon', 'disks', 'k0_sv', ...
    'reson_sv', 'GeomDisks_sv', 'clockout', 'duration',...
    'Rm_sv0','Tm_sv0','Rp_sv0','Tp_sv0','GeomDisks_sv0')

else
 beep; beep; beep
end

if POOL; matlabpool close; end

if COMM
disp('%----------- END: Main_Rows ------------%')
disp('%---------------------------------------%')
end

return

%%%%%%%%%%%%%%%%%%%%
%%% SUBFUNCTIONS %%%
%%%%%%%%%%%%%%%%%%%%

function [Rm,Tm,Rp,Tp,v_vec,k0,x_lim,reson_mkr]...
    = fn_MultiMode_MultiFloe_loop(lam_vec, Param, TankDim, ...
    GeomDisks, Mesh, scatyp, fortyp, COMM)

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
    Param.extra_pts, scatyp, COMM);
   
   %Mesh.r_vec, Mesh.th_vec, Mesh.x_vec, Mesh.y_vec, Mesh.FS_mesh, ...

 end % end loop_lam
 
return
 
function GeomDisks=fn_GeomDisks(Np,TankDim,epsilon,COMM,mxc)

count=1; flg=1;

while and(flg~=0,count<=mxc)
 %%% Location and geometry of disks: [x_c y_c Rad thick]
 for loop=1:Np
  epsx = 2*(rand-0.5)*epsilon; epsy = 2*(rand-0.5)*epsilon;
%   GeomDisks(loop,:) = ...
%       [epsx+TankDim(1)/2 epsy+loop*(TankDim(2)/(Np+1)) 0.495 33e-3];
  GeomDisks(loop,:) = ...
      [epsx+TankDim(1)/2 epsy+(TankDim(2)/2/Np)+(TankDim(2)/Np)*(loop-1) 0.495 33e-3];
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
        
