%% citeph_get_calib_prams.m
%% Author: Timothy Williams
%% Date:   20130724, 10:14:53 CEST
function [dirname,T_target,H_target,type,expt_name] = citeph_get_c39_prams(test_num)

DO_TEST  = 0;
if nargin==0
   test_num = 13;
   DO_TEST  = 1;
end

c39_prams = conc39_testspecs();

dirname     = c39_prams(test_num).dirname;
expt_name   = c39_prams(test_num).name;
T_target    = c39_prams(test_num).period;
H_target    = c39_prams(test_num).wave_height;
type	    = c39_prams(test_num).type;

if strcmp(type,'Regular')
   dirname  = ['regular/' dirname];
else
   dirname  = ['irregular/' dirname];
end

if DO_TEST
   dirname,T_target,H_target,type,expt_name
end
