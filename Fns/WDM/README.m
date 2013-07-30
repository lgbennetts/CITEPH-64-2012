% README gives instructions on preparing data to run WDM.
%
% Wave staff data must be in an array called "data(nsamples,nstaffs)"
% nsamples may be any length > 4095.
% nstaffs = np, must be 3 or larger.
% The data array is called from a file named "s870nn", where nn is the run number, 
% which is set in the overall loop of the main program: WDM.m.
% Change the file names in WDM to suit yourself.
%
% The array (of wave staffs) polar coordinates (R,A) must be set in WDM and
% the maximum wavenumber set at the top of WAVENUMS.m. Usually this is set
% to the Nyquist wavenumber; i.e. 2*pi/(2*diameter of the array of staffs).
%
% Open 2 folders "fspect" and "kspect" to accept the spectra from
% "WAVEPLOTS" and "WAVENUMS" respectively. They are called from WDM.m 
% run by run.
%
% Finally set the frequency limits in WDM: lf and hf and the sampling frequency of the 
% data, ns. Usually you want to set ns to twice the frequency corresponding
% to the Nyquist wavenumber or a little higher. You can also select the
% number of voices, nv, into which each octave is broken. Usually nv = 4.