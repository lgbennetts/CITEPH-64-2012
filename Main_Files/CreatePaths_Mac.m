function CreatePaths_Mac

path(pathdef)

path(path,'..')
path(path,'../Fns')
path(path,'../Fns/Parametrisation')
path(path,'../Fns/AttnModels')
path(path,'../Fns/AttnModels/Roots')
path(path,'../Fns/Misc')
path(path,'../Fns/CollectData')
path(path,'../Fns/CITEPH_specifics')
path(path,'../Fns/SpectraPrelim')
path(path,'../Fns/SpectraPrelim/BumpMethod')
path(path,'../Fns/SpectraPrelim/LookAtWaveSpectra')
path(path,'../Fns/SpectraPrelim/LookAtWaveSpectra/Current')
path(path,'../Fns/Paper')
path(path,'../Fns/Paper/Clean')
path(path,'../Fns/Paper/Collisions')
path(path,'../Fns/Paper/IrregularTransmission')
path(path,'../Fns/SpectraPrelim/OverwashPrelim')
path(path,'../Fns/SkeneOverwash/FloatingPlateSimulation')
path(path,'../Fns/SkeneOverwash/EnergyAnalysis')
path(path,'../Fns/SkeneOverwash/Validation')

path(path,'../Fns/MovingFourier')
path(path,'Temp_data/');

%path(path,'../Fns/Misc/Weight_Norm')
%path(path,'../Fns/Misc/Weight_Ones')
path(path,'../Fns/Misc/Weight_Standard')

if or(strcmp(getenv('LOGNAME'),'a1612881'),...
  strcmp(getenv('LOGNAME'),'lbennetts'))
 %path(path,'../../EXTRA_MATLAB_Fns'); % For colour
 path(path,'../../EXTRA_MATLAB_Fns/GEN_progs');
 %path(path,'../../EXTRA_MATLAB_Fns/l1magic');
 %path(path,'../../EXTRA_MATLAB_Fns/l1magic/Optimization');
 %path(path,'../../EXTRA_MATLAB_Fns/cvx');
 path(path,'Temp_data/');
 %path(path,'../../Fabien_Thesis/MovingFourier/');
 path(path,'../Fns/WDM')
 path(path,'../Fns/MLM')
 path(path,'../Fns/AttnModels/GrafInteraction')
 path(path,'../Fns/Greens_fns')
 path(path,'../Fns/MovingFourier')
 path(path,'../Fns/Papers')
 path(path,'../Fns/Papers/CITEPH-Reg')
end

return
