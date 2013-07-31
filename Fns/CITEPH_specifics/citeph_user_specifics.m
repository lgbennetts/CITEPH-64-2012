function basedir = citeph_user_specifics()
%% citeph_user_specifics.m
%% Author: Timothy Williams
%% Date:   20130731, 09:46:06 CEST
%% basedir contains results_preliminary & calibration_waves directories with data in them

user_name   = getenv('LOGNAME');
if strcmp(user_name,'timill')
   basedir  = '/work/timill/CITEPH-data';
elseif strcmp(user_name,'a1612881')
   basedir  = ['/Volumes/My_Passport/CITEPH_2012-053_CITEPH' ...
    '_WAVE_PROPOGATION_IN_ICE-COVERED_SEAS'];
end
