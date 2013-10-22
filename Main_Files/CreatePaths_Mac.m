function CreatePaths_Mac

%pathroot = '/Users/a1612881/Documents/MatLab/';

path(pathdef)

path(path,'..')
path(path,'../Fns')
path(path,'../Fns/Greens_fns')
path(path,'../Fns/Misc')
path(path,'../Fns/Parametrisation')
path(path,'../Fns/CollectData')
path(path,'../Fns/AttnModels')
path(path,'../Fns/AttnModels/Roots')
path(path,'../Fns/AttnModels/GrafInteraction')
path(path,'../Fns/WDM')
path(path,'../Fns/MLM')
path(path,'../Fns/CITEPH_specifics')
path(path,'../Fns/MovingFourier')

%path(path,'../Fns/Misc/Weight_Norm')
%path(path,'../Fns/Misc/Weight_Ones')
path(path,'../Fns/Misc/Weight_Standard')

if strcmp(getenv('LOGNAME'),'a1612881')
 path(path,'../../EXTRA_MATLAB_Fns'); % For colour
 path(path,'../../EXTRA_MATLAB_Fns/GEN_progs');
 path(path,'Temp_data/');
 path(path,'../../Fabien_Thesis/MovingFourier/');
end

return
