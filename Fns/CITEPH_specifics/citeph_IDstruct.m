%% citeph_data_in.m
%% Author: Timothy Williams
%% Date:   20130715, 10:20:20 CEST
%% Interface for the experiments - sensor locations and matrices with data;

function [id_lhs,id_rhs,freqs,nfft_lhs,nfft_rhs,wave_TH,Hs_ind] = citeph_IDstruct()

GET_DATA = 1;
DO_PLOT  = 0;

%%tank (002.pdf)
tank_width  = 2000+4000*3+2000;%mm
tank_len    = 24000;%mm - wavemaker to beach at end??
tank_dep    = 3000;%mm

%%xlim and ylim of tank;
ab = tank_len*[-.5 .5];
cd = tank_width*[-.5 .5];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%(x,y) coords of sensors:
[xy_lhs,xy_rhs]   = citeph_sensor_spots();%m

%%z position (at surface)
x_lhs    = 1e3*xy_lhs(:,2);%mm
xyz_lhs  = [xy_lhs,tank_dep/1e3+0*x_lhs];%m
x_rhs    = 1e3*xy_rhs(:,2);%mm
xyz_rhs  = [xy_rhs,tank_dep/1e3+0*x_rhs];%m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%**************************
%%***********LHS************
%%**************************


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%sensor types:
dtype_lhs  = {'elev' 'elev'   'elev'  'elev'  'elev'  'elev' 'elev'   'elev'  'elev'};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% read in data here:
if GET_DATA==0
   data_lhs          = [];
   sample_rate_lhs   = [];
else
   test_num                         = 1
   [data_lhs,sample_rate_lhs,freqs,wave_TH,Hs_ind0]  =...
         citeph_testspec(test_num,xyz_lhs);
   nfft_lhs                         = size(data_lhs,1);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if 1
   j_pent      = 3:8;
   xyz_lhs     = xyz_lhs(j_pent,:);
   data_lhs    = data_lhs(:,j_pent);
   dtype_lhs   = dtype_lhs(j_pent);
   ONLY_PENT   = 1;
else
   ONLY_PENT   = 0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
id_lhs   = struct('data',data_lhs,...
                  'layout',xyz_lhs',...
                  'datatypes',{dtype_lhs},...
                  'depth',tank_dep/1e3,...
                  'fs',sample_rate_lhs);
%id_lhs   = struct('datatypes',{dtype_lhs})



%size(id_lhs.data)
%id_lhs.layout
%id_lhs.datatypes{9}
%id_lhs.depth
%id_lhs.fs
%return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%**************************
%%***********RHS************
%%**************************

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%sensor types:
dtype_rhs  = {'elev' 'elev' 'elev'   'elev'  'elev'  'elev'  'elev' 'elev'   'elev'  'elev'};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% read in data here:
if GET_DATA==0
   data_rhs          = [];
   sample_rate_rhs   = [];
else
   [data_rhs,sample_rate_rhs,freqs,wave_TH,Hs_ind] =...
         citeph_testspec(test_num,xyz_rhs);
   nfft_rhs                         = size(data_rhs,1);
   Hs_ind   = [Hs_ind0;Hs_ind];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ONLY_PENT
   j_pent      = 3:8;
   Hs_ind      = Hs_ind(:,j_pent);
   xyz_rhs     = xyz_rhs(j_pent,:);
   data_rhs    = data_rhs(:,j_pent);
   dtype_rhs   = dtype_rhs(j_pent);
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
id_rhs   = struct('data',data_rhs,...
                  'layout',xyz_rhs',...
                  'datatypes',{dtype_rhs},...
                  'depth',tank_dep/1e3,...
                  'fs',sample_rate_rhs);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if DO_PLOT
   figure(1);
   clf;

   lab_lhs  = {'S1' 'S2'   'S3'  'S4'  'S5'  'S6' 'S7'   'S8'  'S9'};
   lab_rhs  = {'S10' 'S11' 'S12' 'S13' 'S14' 'S15' 'S16' 'S17' 'S18'};

   for j=1:9
      hold on;
      x  = xyz_lhs(j,1)*1e3;
      y  = xyz_lhs(j,2)*1e3;
      plot(x,y,'x');
      txt   = text(x,y,lab_lhs{j});
      %%
      x  = xyz_rhs(j,1)*1e3;
      y  = xyz_rhs(j,2)*1e3;
      plot(x,y,'x');
      txt   = text(x,y,lab_rhs{j});
      hold off;
   end
   xlim(ab);
   ylim(.2*cd);
end
