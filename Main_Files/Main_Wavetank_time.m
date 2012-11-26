function Main_Wavetank_time

% Main_Wavetank_time
%
% This is the main program used to solve the transient response of a set of
% a set of floatings disks in a 3D wavetank.

disp('%---------- START: Main_Wavetank_time -----------%')

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

%%% Time to frequency domain transformation
[Signal,Sample] = Time2Freq(Param.g,TankDim(3));

%%% Generation of horizontal mesh points for the free-surface region
Mesh = Mesh_FS_def(GeomDisks, TankDim);

%%% Modal parameters and accuracy
Param = ModParam_def(Param);

%% Frequency domain response
addpath('Frequency Solver')
Res_mat_freq = Freq_response_gen_multi(GeomDisks,TankDim,Sample,Param);

%% Generation of results
addpath('Results routines')
Results_time(Res_mat_freq,Signal,Sample,Param,Plate,Tank);


disp('%---------- END: Main_Wavetank_time -----------%')