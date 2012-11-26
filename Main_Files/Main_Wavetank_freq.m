function Main_Wavetank_freq

% Main_Wavetank_freq
%
% This is the main program used to solve the frequency response of a set of
% floating disks to wave maker forcing in a 3D wave tank.

disp('%---------- START: Main_Wavetank_freq -----------%')

%% Creation of local paths
CreatePaths_Mac

%% Parametrisation
% Note that the origin of the Cartesian system is located at the bottom 
% left corner of the wave tank and coincides with the free-surface at rest.

%%% Location and geometry of disks: [x_c y_c Rad thick]
GeomDisks=[4 2.5 .1 1e-4];

%%% Geometry of the wave tank: [length width height]
TankDim = [10 5 2]; 

%%% Check disk/disk and disk/wall overlaps
overlap(GeomDisks, TankDim)

%%% Physical properties
Param = ParamDef3d(GeomDisks);

%%% Forcing
Forcing = Force_def(Param.g, TankDim(3));

%%% Generation of horizontal mesh points for the free-surface region
Mesh = Mesh_FS_def(GeomDisks, TankDim);

%%% Modal parameters and accuracy
Param = ModParam_def(Param);

%% Solution - Disks in a channel
% Scattering matrix is S = [Rm Tp;Tm Rp]
% Wavenumbers matrices v_vec and u_vec are scaled by vertical wavenumbers.
[Rm,Tm,Rp,Tp,v_vec,u_vec,k0,weight_0,x_lim] = ...
    fn_MultiMode_MultiFloe([Param.rho_0 Param.rho Param.nu Param.E ...
    Param.rho_0], Param.N, Param.Mev, TankDim, Forcing.kappa, ...
    Param.beta.', GeomDisks(:,4).', Param.draft.', GeomDisks(:,3).', ...
    [GeomDisks(:,1).';GeomDisks(:,2).'], Mesh.r_vec, Mesh.th_vec, ...
    Mesh.x_vec, Mesh.y_vec, Mesh.FS_mesh, Param.res_green, ...
    Param.extra_pts);


%% Solution - Wave tank 
[Amp_1P,Amp_1M,Amp_3P,Amp_3M] = ...
    Wavetank3D_solver(Forcing.AmpWM, k0, weight_0, v_vec, Rm, Tm, Rp, ...
    Tp, TankDim(2),TankDim(3), [0 x_lim TankDim(1)], 2*pi*Forcing.f, ...
    Param.N-1, size(u_vec,2)-1, size(u_vec,2) - Param.Mev, Param.Rb, 1);

disp('%----------- END: Main_Wavetank_freq ------------%')

return