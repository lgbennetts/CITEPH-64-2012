function Param = ModParam_def(Param,extra_pts)

% Param = ModParam_def(Param)
%
% Setup of modal parameters and accuracy parameters for the solution
% method.

%%% Number of vertical modes (travel + evan)
Param.N = 1;

%%% Number of evanescent horizontal modes
Param.Mev = 0;

%%% Number of extra points for irregular frequencies
if ~exist('extra_pts','var')
 Param.extra_pts = 0;
else
 if isempty(extra_pts)
  Param.extra_pts = 0;
 else
  Param.extra_pts = extra_pts;
 end
end

%%% Number of points for the intergration of the Green's function
Param.res_green = 100;

%%% Max number of terms used in Green's function
Param.terms_green = 100;

%%% Limit for cutoff of Greens fn series
Param.cutoff_green = 1e-1;

%%% Tolerance on resonances

Param.tolres=1e-5; %7.5e-4; %0; %

return