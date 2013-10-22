function basedir = citeph_user_specifics()
%% citeph_user_specifics.m
%% Author: Timothy Williams
%% Date:   20130731, 09:46:06 CEST
%%
%% basedir contains the directories
%% 'results_preliminary' & 'calibration_waves' directories with raw data in them

[user_number,user_name]  = citeph_get_user;

switch user_number
case 1%%Tim
   [comp_number,comp_name] = TW_which_computer;
   switch comp_number
   case 1%%desktop
      basedir  = '/Volumes/Tim_Ext_HD2/WORK/Data/EXPERIMENTS/CITEPH/RAW_data';
   case 2%%laptop
      basedir  = '/Volumes/Tim_Ext_HD2/WORK/Data/EXPERIMENTS/CITEPH/RAW_data';
   case 3%%hexagon
      basedir  = '/work/timill/CITEPH-data/RAW_data';
   case 4%%nansen
      basedir  = ''
   end

case 2%%Luke
   basedir  = ['/Volumes/My_Passport/CITEPH_2012-053_CITEPH' ...
    '_WAVE_PROPOGATION_IN_ICE-COVERED_SEAS'];
end

if strcmp(basedir,'')
   disp('');
   disp('Warning (citeph_user_specifics.m): no directory <<basedir>> defined')
   disp('');
end
