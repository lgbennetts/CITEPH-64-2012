function Main_Channel_freq(Np, lam_vec,extra_pts,file_marker)

if ~exist('Np','var'); Np=1; end
if ~exist('extra_pts','var'); extra_pts=[]; end
if ~exist('file_marker','var'); file_marker=0; end
if ~exist('lam_vec','var'); lam_vec=1.0*pi/16; end

% Main_Wavetank_freq
%
% This is the main program used to solve the frequency response of a set of
% floating disks to wave maker forcing in a 3D wave tank.

disp('%---------- START: Main_Channel_freq -----------%')

%% Creation of local paths
%CreatePaths_Mac

%% Wavelengths
%%% NB. There are m `tank modes' when k\in \pi(m,m+1)/width

%lam_vec = 1:0.025:3;

%% Parametrisation
% Note that the origin of the Cartesian system is located at the bottom 
% left corner of the wave tank and coincides with the free-surface at rest.

%%% Geometry of the wave tank: [length width height]
TankDim = [10 16 2]; 

%%% Location and geometry of disks: [x_c y_c Rad thick]

GeomDisks = zeros(Np, 4);

for loop=1:Np
 GeomDisks(loop,:)=[4 0+(TankDim(2)/2/Np)+(TankDim(2)/Np)*(loop-1) 0.5 1e-10];
end

%%% Check disk/disk and disk/wall overlaps
overlap(GeomDisks, TankDim)

%%% Physical properties
Param = ParamDef3d(GeomDisks);

%%% Generation of horizontal mesh points for the free-surface region
Mesh = Mesh_FS_def(GeomDisks, TankDim);

%%% Modal parameters and accuracy
Param = ModParam_def(Param,extra_pts);

for loop_lam=1:length(lam_vec)

%%% Forcing
Forcing = Force_def(Param.g(1), TankDim(3), 'waveno', lam_vec(loop_lam));

%% Solution - Disks in a channel
% Scattering matrix is S = [Rm Tp;Tm Rp]
% Wavenumbers matrices v_vec and u_vec are scaled by vertical wavenumbers.
[Rm,Tm,Rp,Tp,v_vec,u_vec,k0,weight_0,x_lim] = ...
    fn_MultiMode_MultiFloe([Param.rho_0 Param.rho Param.nu Param.E ...
    Param.g], Param.N, Param.Mev, TankDim, Forcing.kappa, ...
    Param.beta.', GeomDisks(:,4).', Param.draft.', GeomDisks(:,3).', ...
    [GeomDisks(:,1).';GeomDisks(:,2).'], ...
    Mesh.r_vec, Mesh.th_vec, Mesh.x_vec, Mesh.y_vec, Mesh.FS_mesh, ...
    [Param.res_green, Param.terms_green, Param.cutoff_green, Param.tolres], ...
    Param.extra_pts);

R_vec{loop_lam} = Rm; T_vec{loop_lam} = Tm; 

v_vecs{loop_lam} = v_vec(1,:); %k_vec(loop_lam) = k0(1);

end

if and(~file_marker, length(lam_vec)==1)
 disp('Rm='); disp(abs(Rm))
 disp('Tm='); disp(abs(Tm))
%  disp('Rp='); disp(abs(Rm))
%  disp('Tp='); disp(abs(Tm))
end

if file_marker~=0
    
 w = TankDim(2);
 
 for loopi=1:4
  dum_str{loopi} = num2str(GeomDisks(1,loopi));    
  for loopii=2:Np   
    dum_str{loopi} = [ dum_str{loopi}, ', ' ...
        num2str(GeomDisks(loopii,loopi)) ];
  end
 end
 clear loopi loopii
 dum_str{5} = num2str(Param.E(1));    
  for loopii=2:Np   
    dum_str{5} = [ dum_str{5}, ', ' ...
        num2str(Param.E(loopii)) ];
  end
  clear loopii 

 description = {['Tank: ' num2str(TankDim(1)) 'x' num2str(TankDim(2)) ...
     'x' num2str(TankDim(3))], ...
     [int2str(Np) ' Disk(s): x=' dum_str{1} ...
     '; y=' dum_str{2} ...
     '; r=' dum_str{3} ...
     '; h=' dum_str{4} ...
     '; E=' dum_str{5} ], ...
     ['Vert modes=' int2str(Param.N) ...
     ', extra points=' int2str(Param.extra_pts) ...
     ', angular res=' int2str(Param.res_green) ...
     ', Green fn terms=' int2str(Param.terms_green) ...
     ', Greens fn cutoff=' int2str(Param.cutoff_green) ] };
     
 file_name = ['Main_Channel_freq_',file_marker,'.mat'];

 save(['../../../../../Documents/MatLab/Data/3d_Wavetank/Tests/',file_name],...
    'R_vec', 'T_vec', 'lam_vec', 'v_vecs', 'w', 'description')

else
 beep; beep; beep
end

disp('%----------- END: Main_Channel_freq ------------%')

return