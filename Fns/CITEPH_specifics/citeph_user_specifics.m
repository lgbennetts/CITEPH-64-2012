function basedir = citeph_user_specifics()
%% citeph_user_specifics.m
%% Author: Timothy Williams
%% Date:   20130731, 09:46:06 CEST
%% basedir contains results_preliminary & calibration_waves directories with data in them

user_name   = getenv('LOGNAME');
if strcmp(user_name,'timill')
   basedir  = '/work/timill/CITEPH-data/results_preliminary/';
   cprintf('blue','Tim: your basedir needs modification\n'); 
elseif strcmp(user_name,'a1612881')
   basedir  = '/Volumes/Scratch/Data/CITEPH-64-2012/';
end

return