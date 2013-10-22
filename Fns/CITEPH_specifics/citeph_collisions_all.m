function citeph_collisions_all

%% in citeph_coll0 below;
%% adjust DO_CALC & DO_ANALYSIS to speed up or redo calculations;
%% adjust SENSOR_LOOP to change which direction(s) of acceleration you want

tt = {'c39','c79'};
for j=2
   test_type   = tt{j};
   citeph_coll0(test_type)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function citeph_coll0(test_type)

tres  = .004;%%time resolution [s]

if strcmp(test_type,'c79')
   Ntests      = 19;
   Ntests_reg  = 15;
elseif strcmp(test_type,'c39')
   Ntests      = 12;
   %Ntests_reg  = 8;
elseif strcmp(test_type,'calib')
   Ntests   = 17;
end

DO_CALC     = 1;%%if need to do or redo calculation set to 1; else if it's done set to 0;
DO_ANALYSIS = 1;%%if it's done (& the files saved) can speed things up by setting this to 0;
if DO_CALC==0
   figure(14);
   rms_a    = zeros(Ntests_reg,1);
   T_coll   = zeros(Ntests_reg,1);
   Hi_Frac  = zeros(Ntests_reg,1);
end

%%***********************
%%******SENSOR LOOP******
%%***********************
opts     = {'Ax','Ay','Az'};
%for n=1:3
for n=1
   opt      = opts{n};
   outfile  = ['out/coll_stats_' opt '.mat'];
   figfile  = ['out/collision_strength_' opt '.eps'];
   %%
   if DO_ANALYSIS | DO_CALC
      for test_num=1:Ntests_reg
      %for test_num=1:Ntests
      %for test_num=6
         if DO_CALC%%do calculation
            citeph_coll_1test(test_type,test_num,opt)
         else%%open results & plot;
            [Tcol,Acol,Asig,Tp,Hs,sensor_names,Col_dt] =...
               citeph_coll_1test_OpenResultsFile(test_type,test_num,opt);
            T_target(test_num,1) = Tp;
            H_target(test_num,1) = Hs;

            Nsens = length(Tcol);
            for r=1:Nsens
               acol2 = Acol{r};
               
               %%collision period (mean length of time between collisions);
               tcol2                = Tcol{r};
               dt                   = tcol2(2:end)-tcol2(1:end-1);
               T_coll(test_num,r)   = mean(dt);%%6x1 vector
               dt_all{test_num,r}   = dt;
               %%
               col_dt(test_num,r)   = mean(Col_dt{r});
               %%
               rms_a(test_num,r)    = sqrt(mean(acol2.^2));%%6x1 vector

               %%fraction of variance that is high frequency (percentage)
               asig1                = Asig(r,1);
               asig2                = Asig(r,2);
               Hi_Frac(test_num,r)  = (asig2/asig1)^2*100;%%6x1 vector
               %pause
            end
            if ~exist('out')
               eval(['!mkdir out']);
            end
            save(outfile,'rms_a','T_target','H_target','sensor_names','T_coll','Hi_Frac','dt_all','col_dt');
            %plot(T_target,rms/H_target,'x');
            %hold on
         end
      end
   else%%do/don't do calc
      load(outfile);
      %T_coll
      %dt_all
   end
end%n - opt (Ax,Ay,Az)

if DO_CALC==0
   Nsens = size(rms_a,2);
   if strcmp('Az',opt)
      colin = {'xb','xr','or','xk'};
   else
      colin = {'xb','ob','ok','xr','or','xk'};
   end
   %%
   subplot(3,1,1);
   for n=1:Nsens
      P(n)  = plot(T_target,rms_a(:,n)./H_target,colin{n});
      hold on;
   end
   legend(P,sensor_names,'location','eastoutside');
   hold off;
   GEN_proc_fig('T_p, s','rms(a_{col})/H_s, s^{-2}');
   %%
   subplot(3,1,2);
   for n=1:Nsens
      P(n)  = plot(T_target,T_coll(:,n)./T_target,colin{n});
      hold on;
   end
   legend(P,sensor_names,'location','eastoutside');
   ylim([0 3.5]);
   hold off;
   GEN_proc_fig('T_p, s','T_{col}/T_p');
   %%
   subplot(3,1,3);
   for n=1:Nsens
      P(n)  = plot(T_target,Hi_Frac(:,n),colin{n});
      hold on;
   end
   legend(P,sensor_names,'location','eastoutside');
   hold off;
   GEN_proc_fig('T_p, s','High-frequency fraction of variance, %');
   %%
   saveas(gcf,figfile,'epsc');
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   loc   = [5 6 1 3 4 2];
   %for r=1:Nsens
   for r=4
      N1          = size(dt_all,1);
      Dt_all      = [];
      check_ind   = 1;
      for test_num=1:N1
         Tp       = T_target(test_num);
         Hs       = H_target(test_num);
         %Dt_all   = [Dt_all;dt_all{test_num,r}/Tp];
         if check_ind
            Dt_all   = dt_all{test_num,r};
            Y        = Dt_all/Tp;
            Yav      = mean(Y);
            figure(15);
            subplot(3,2,loc(r));
            hist(Y,50);
            ttl   = title([sensor_names{r},': ',num2str(Tp),' s, ',num2str(Hs),'m; mean = ',num2str(Yav,'%0.1f')]);
            GEN_font(ttl);
            hold on;
            yl = get(gca,'ylim');
            plot(Yav+0*yl,yl,'g');
            ylim(yl);
            hold off;

            GEN_proc_fig('\Delta t/T_p','frequency');
            pause
         else%% bunch them together
            Dt_all   = [Dt_all;dt_all{test_num,r}];
         end
      end
      % {r,loc(r),sensor_names{r}}
      if ~check_ind
         figure(15);
         subplot(3,2,loc(r));
         hist(Dt_all,50);
         ttl   = title(sensor_names{r});
         GEN_font(ttl);
         GEN_proc_fig('\Delta t/T_p','frequency');
      end
   end
end


function citeph_coll_1test(test_type,test_num,opt)

if nargin==0
   test_type   = 'c79';
   test_num    = 19;
end

[time,data,file_list,sensor_names]  = citeph_get_data(test_type,test_num,opt,'true');
nsens = 6;

if 0
   figure(1);data(:,nsens)
   plot(time(:,6),data(:,nsens))
   pause
end

if strcmp(test_type,'c79')
   [expt_dir,T_target,H_target,type,expt_name] = citeph_get_c79_prams(test_num,'true');
   tt2   = 'conc_79/';
else%% 'c39'
   [expt_dir,T_target,H_target,type,expt_name] = citeph_get_c39_prams(test_num,'true');
   tt2   = 'conc_39/';
end

Nsens = size(data,2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
basedir  = citeph_user_specifics;
jj       = find(basedir=='/');
basedir2 = [basedir,'/../PROCESSED_data/'];
if ~exist(basedir2)
   eval(['!mkdir ' basedir2]);
end
basedir2 = [basedir2,tt2];
if ~exist(basedir2)
   eval(['!mkdir ' basedir2]);
end

if strcmp(type,'Regular')
   type  = 'regular';
else
   type  = 'irregular';
end
basedir3 = [basedir2,type,'/'];%%/regular or /irregular
if ~exist(basedir3)
   eval(['!mkdir ' basedir3]);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for n=1:Nsens
%for n=nsens
   time0 = time(:,n);
   displ = data(:,n);
   inloc = file_list{n}
   expt_dir
   %%
   outloc                  = {basedir2,expt_dir,opt,sensor_names{n} }
   %[an_mat,Sn_mat,t_int]   = citeph_1sensor_collisions(time0,displ,T_target,outloc,inloc);
   citeph_1sensor_collisions(time0,displ,T_target,H_target,outloc,inloc);
end

function [Tcol,Acol,Asig,T_target,H_target,sensor_names,Col_dur] =...
   citeph_coll_1test_OpenResultsFile(test_type,test_num,opt)

if nargin==0
   test_type   = 'c79';
   test_num    = 19;
end

[time,data,file_list,sensor_names]  = citeph_get_data(test_type,test_num,opt,'true');
nsens = 6;

if 0
   figure(1);data(:,nsens)
   plot(time(:,6),data(:,nsens))
   pause
end

if strcmp(test_type,'c79')
   [expt_dir,T_target,H_target,type,expt_name] = citeph_get_c79_prams(test_num,'true');
   tt2   = 'conc_79/';
else%% 'c39'
   [expt_dir,T_target,H_target,type,expt_name] = citeph_get_c39_prams(test_num,'true');
   tt2   = 'conc_39/';
end

Nsens = size(data,2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
basedir  = citeph_user_specifics;
jj       = find(basedir=='/');
basedir2 = [basedir,'/../PROCESSED_data/'];
basedir3 = [basedir2,tt2,expt_dir,'/',opt,'/mat_files/'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for n=1:Nsens
%for n=nsens
   time0 = time(:,n);
   displ = data(:,n);
   inloc = file_list{n};
   fname = [basedir3,'/collisions_',sensor_names{n},'.mat']
   load(fname);
   Tcol{n}     = tcol2;
   Acol{n}     = acol2;
   Asig(n,:)   = [asig1,asig2];
   Col_dur{n}  = col_dt;
end
