%% citeph_get_calib_prams.m
%% Author: Timothy Williams
%% Date:   20130724, 10:14:53 CEST
function [dirname,T_target,H_target,type,expt_name] = citeph_get_c79_prams(test_num,scale)

DO_TEST  = 0;
if nargin==0
   test_num = 19;
   DO_TEST  = 1;
elseif 0
   test_num = 13;
   DO_TEST  = 1;
end

if ~exist('scale')
   scale = 'false';
end
c79_prams   = conc79_testspecs();

dirname     = c79_prams(test_num).dirname;
expt_name   = c79_prams(test_num).name;
T_target    = c79_prams(test_num).period;
H_target    = c79_prams(test_num).wave_height;
type	    = c79_prams(test_num).type;

if strcmp(type,'Regular')
   dirname  = ['regular/' dirname];
else
   dirname  = ['irregular/' dirname];
end

if strcmp(scale,'true')%%scale variables to measured/basin scale;
   T_target = T_target/10;
   H_target = H_target/100;
end

if DO_TEST
   dirname,T_target,H_target,type,expt_name
end
