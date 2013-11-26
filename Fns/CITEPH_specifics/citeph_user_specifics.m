function basedir = citeph_user_specifics(tests)
%% citeph_user_specifics.m
%% Author: Timothy Williams
%% Date:   20130731, 09:46:06 CEST
%%
%% basedir contains the directories
%% 'results_preliminary' & 'calibration_waves' directories with raw data in them

if ~exist('tests','var'); tests='final'; end

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
  dum = {'/Volumes/Scratch/Data/','/Volumes/My_Passport/'};
  eval(['ex=exist(dum{1});'])
  count=1;
  while ex~=7
   count=count+1;
   if count<length(dum)
    cprintf('green',['>>> ' dum ' not found\n'])
    s = input('Enter new base dir or <return> to exit', 's');
    if isempty(s)
     return
    else
     eval(['ex=exist(s);'])
    end
   else
    eval(['ex=exist(dum{count});'])
   end
  end
  if strcmp(tests,'prelim')
   basedir  = [dum{count} 'CITEPH-64-2012/CITEPH-data-preliminary/'];
  elseif strcmp(tests,'final')
   basedir  = [dum{count} 'CITEPH-64-2012/CITEPH-data-final/'];
  end
end

if strcmp(basedir,'')
 disp('');
 disp('Warning (citeph_user_specifics.m): no directory <<basedir>> defined')
 disp('');
end
