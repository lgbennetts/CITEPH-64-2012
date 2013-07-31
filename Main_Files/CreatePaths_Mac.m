function CreatePaths_Mac

%pathroot = '/Users/a1612881/Documents/MatLab/';

path(pathdef)

path(path,'..')
path(path,'../Fns')
path(path,'../Fns/Greens_fns')
path(path,'../Fns/Misc')
path(path,'../Fns/Parametrisation')
path(path,'../Fns/CollectData')
path(path,'../Fns/Floe')
path(path,'../Fns/WDM')
path(path,'../Fns/MLM')
path(path,'../Fns/CITEPH_specifics')
path(path,'../Fns/MovingFourier')

path(path,'../Fns/Misc/Weight_Norm')

if strcmp(getenv('LOGNAME'),'a1612881')
 path(path,'../../EXTRA_MATLAB_Fns'); % For colour
 path(path,'../../EXTRA_MATLAB_Fns/GEN_progs');
end

return
