function [xy_lhs,xy_rhs] = citeph_sensor_spots()

%% citeph_data_in.m
%% Author: Timothy Williams
%% Date:   20130715, 10:20:20 CEST
%% Interface for the experiments - sensor locations and matrices with data;

%%%tank (002.pdf)
%tank_width  = 2000+4000*3+2000;%mm
%tank_len    = 24000;%mm - wavemaker to beach at end??
%tank_dep    = 3000;%mm
%
%%%xlim and ylim of tank;
%ab = tank_len*[-.5 .5]
%cd = tank_width*[-.5 .5]


%%**************************
%%***********LHS************
%%**************************

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% locations of sensors on lhs
%% detail A (001.pdf):
x_pent   = [770   943   1222  1222  943];
xcp   = mean(x_pent);%%centre of pentagon (S19)

%%                      | ------- pentagon ------- |
lab_lhs  = {'S1' 'S2'   'S3'  'S4'  'S5'  'S6' 'S7'   'S8'  'S9'  'S19'};
x_lhs0   = [0     160   770   943   1222  1222  943   1710  2840  1020]';
y_lhs0   = [0     0     0     238   147   -147  -238  0     0     0   ]';
   %% relative to S1

%001.pdf: positions of S1 relative to central axis
x_s1  = -8340;
y_s1  = 0;

x_lhs    = x_s1+x_lhs0;
y_lhs    = y_s1+y_lhs0;
xy_lhs  = [x_lhs y_lhs]/1e3;%m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%**************************
%%***********RHS************
%%**************************

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% locations of sensors on rhs
%% detail B (001.pdf):
x_pent   = [770   943   1222  1222  943];
xcp   = mean(x_pent);%%centre of pentagon (S20)

%%                      | ------- pentagon ------- |
lab_rhs  = {'S10' 'S11' 'S12' 'S13' 'S14' 'S15' 'S16' 'S17' 'S18' 'S20'};
x_rhs0   = [0     160   770   943   1222  1222  943   1710  2840  1020]';
y_rhs0   = [0     0     0     238   147   -147  -238  0     0     0   ]';
   %% relative to S10

%001.pdf: positions of S18 relative to central axis
% (get S10 relative to S18 from detail B)
x_s18 = 8340;
x_s10 = x_s18+(x_rhs0(1)-x_rhs0(end));
y_s10 = 0;


x_rhs    = x_s10+x_rhs0;
y_rhs    = y_s10+y_rhs0;
xy_rhs   = [x_rhs y_rhs]/1e3;%m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
