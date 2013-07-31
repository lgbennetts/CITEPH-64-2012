function citeph_batch_fourier_series(test_type,data_type,opt)
%% citeph_batch_fourier_series.m
%% Author: Timothy Williams
%% Date:   20130724, 13:36:06 CEST
%% test_num    = test no
%% test_type   = 'calib','c79' or 'c39'
%% data_type   = 'S','S_zoom','Ax','Ay','Az'
%% opt         = 1 (moving fft), 2 (look for collisions), 3 (single fft over full record)
if nargin==0
   test_type   = 'c79';
   %data_type   = 'Ax';
   %data_type   = 'Ay';
   data_type   = 'Az';
   %data_type   = 'S';
   opt         = 1;
end
if strcmp(test_type,'c79')
   Ntests   = 19;
else strcmp(test_type,'c39')
   Ntests   = 12;
else
   Ntests   = 17;
end

%for test_num=1:Ntests
for test_num=1
   batch_fourier_series(test_num,test_type,data_type,opt);
end

function batch_fourier_series(test_num,test_type,data_type,opt)
%% test_num    = test no
%% test_type   = 'calib','c79' or 'c39'
%% data_type   = 'S','S_zoom','Ax','Ay','Az'
%% opt         = 1 (moving fft), 2 (look for collisions), 3 (single fft over full record)

if nargin==0
   test_num    = 1;
   test_type   = 'c79';
   %data_type   = 'Ax';
   %data_type   = 'Ay';
   data_type   = 'Az';
   %data_type   = 'S';
   opt         = 1;
end

if strcmp(test_type,'calib')
   [dirname,T_target,H_target,type,expt_name] = citeph_get_calib_prams(test_num);%%NB full scale
elseif strcmp(test_type,'c79')
   [dirname,T_target,H_target,type,expt_name] = citeph_get_c79_prams(test_num)%%NB full scale
else%%'c39'
   [dirname,T_target,H_target,type,expt_name] = citeph_get_c39_prams(test_num);%%NB full scale
end

%[time,data,file_list]   = citeph_get_data(test_type,test_num,data_type);%%NB full scale
[time,data,file_list]   = citeph_get_data_v2(test_type,test_num,data_type);%%NB full scale

if 1%% scale data
   scale_fac   = 100;
   %%
   time        = time/sqrt(scale_fac);
   T_target    = T_target/sqrt(scale_fac);
   %%
   data        = data/scale_fac;
   H_target    = H_target/scale_fac;
end

j_lhs    = [1:10];
j_rhs    = [11:20];
Nprobes  = length(file_list)
%%
outloc      = cell(4,1);
jb          = find(file_list{1}=='/');
outloc{1}   = [file_list{1}(1:jb(end-1)-1),'_processed'];
   %% eg /work/timill/CITEPH-data/results_preliminary/conc_79/regular_processed
outloc{2}   = file_list{1}(jb(end-2)+1:jb(end-1)-1);%% eg 23051525.a13
outloc{3}   = ['/' data_type '/'];

if 0
   file_list{1}
   for j=1:3
      outloc{j}
   end
   return;
end

if strcmp(data_type,'Az')
   sname_list  = {'A1z';'A4z';'A5z';'A6z'};
elseif strcmp(data_type,'Ax')
   sname_list  = {'A1x';'A2x';'A3x';'A4x';'A5x';'A6x'};
elseif strcmp(data_type,'Ay')
   sname_list  = {'A1y';'A2y';'A3y';'A4y';'A5y';'A6y'};
elseif strcmp(data_type,'S')
   for n=1:9
      sname_list{n}  = ['S0' num2str(n)];
   end
   sname_list{10}  = 'S19';
   for n=11:19
      sname_list{n}  = ['S' num2str(n-1)];
   end
   sname_list{20}  = 'S20';
else
   for n=1:9
      sname_list{n}  = ['S0' num2str(n) '_zoom'];
   end
   sname_list{10}  = 'S19_zoom';
   for n=11:19
      sname_list{n}  = ['S' num2str(n-1) '_zoom'];
   end
   sname_list{20}  = 'S20_zoom';
end

for n=1:Nprobes
   sensor_name = sname_list{n};
   outloc{4}   = sensor_name;
   displ       = data(:,n);
   if strcmp(data_type,'S_zoom')
      time0 = time(:,n);
   else
      time0 = time;
   end

   inloc = file_list{n};
   if opt==1%%moving FFT;
      [an_mat,Sn_mat,t_int]   = citeph_1sensor_movingFFT(time0,displ,T_target,outloc,inloc);
   elseif opt==2%%look for collisions;
      [t_col,a_col]  = citeph_1sensor_collisions(time0,displ,T_target,outloc,inloc);
   else%%single FFT over the whole record;
      [an,Sn]   = citeph_1sensor(time0,displ,T_target,outloc,inloc);
   end
end
