function Param = ModParam_def(Param)

% Param = ModParam_def(Param)
%
% Setup of modal parameters and accuracy parameters for the solution
% method.

%%% Number of evanescent vertical modes
Param.N = 3;

%%% Number of evanescent horizontal modes
Param.Mev = 0;

%%% Number of extra points for irregular frequencies
Param.extra_pts = 0;

%%% Number of points for the intergration of the Green's function
Param.res_green = 100;