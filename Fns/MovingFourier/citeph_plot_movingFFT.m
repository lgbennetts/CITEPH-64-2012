%% citeph_plot_movingFFT.m
%% Author: Timothy Williams
%% Date:   20130729, 10:43:04 CEST
function [PP,TM_vec,Sn_mat] = citeph_plot_movingFFT(test_type,test_num,sensor_name)

if nargin==0
   test_type   = 'c79';
   test_num    = 1;
   %sensor_name = 'S19';
   sensor_name = 'A1z';
end

basedir  = citeph_user_specifics;
if strcmp(test_type,'calib')
   dir0  = [basedir '/calibration_waves/'];
   [dirname,T_target,H_target,type,expt_name] = citeph_get_calib_prams(test_num);%%NB full scale
elseif strcmp(test_type,'c79')
   dir0  = [basedir '/results_preliminary/conc_79/'];
   [dirname,T_target,H_target,type,expt_name] = citeph_get_c79_prams(test_num)%%NB full scale
else%%'c39'
   dir0  = [basedir '/results_preliminary/conc_39/'];
   [dirname,T_target,H_target,type,expt_name] = citeph_get_c39_prams(test_num);%%NB full scale
end

if strcmp(sensor_name(1),'A')
   dir1        = [dir0 '/' dirname '_processed/' sensor_name([1 3])];
   data_type   = 2;
elseif strcmp(sensor_name(1),'S')
   dir1     = [dir0 '/' dirname '_processed/' 'S'];
   data_type   = 1;
elseif strcmp(sensor_name(end-3:end),'zoom')
   dir1     = [dir0 '/' dirname '_processed/' 'S_zoom'];
   data_type   = 1;
end
matfile  = [dir1 '/' 'movingFFT_' sensor_name '.mat']
fig_dir  = [dir1 '/eps_files'];
if ~exist(fig_dir)
   eval(['!mkdir ' fig_dir]);
end
figname1 = [fig_dir  '/movingFFT_' sensor_name '_TpHs.eps'];
figname2 = [fig_dir  '/movingFFT_' sensor_name '.eps'];

load(matfile);

Nwindows = size(an_mat,2);
Nw       = size(an_mat,1);
tm_vec   = mean(t_int,2);
%%
scale_fac   = 100;
T_target    = T_target/sqrt(scale_fac);
H_target    = H_target/scale_fac;

A           = load(inloc);
time        = A(:,1)/sqrt(scale_fac);
displ       = A(:,2)/scale_fac;

DO_PLOT  = 1;
if DO_PLOT
   %%plot time series:
   figure(101);
   subplot(2,1,1);
   plot(tm_vec,Tp_vec,'k');
   hold on;
   plot(tm_vec,Tm02_vec,'--r');
   hold off;
   GEN_proc_fig('time, s','T_p & T_m02, s');
   jb    = find(dirname=='/');
   tstr  = [dirname(jb+1:end),', ' sensor_name];
   ttl   = title(tstr);
   GEN_font(ttl);

   subplot(2,1,2)
   if data_type==2
      Hs_vec   = Hs_vec/2;
      ylab     = 'a_s, ms^{-2}';
   else
      ylab     = 'H_s, m';
   end
   plot(tm_vec,Hs_vec);
   GEN_proc_fig('time, s',ylab);
   saveas(gcf,figname1,'epsc');
   disp(['saved: ' figname1]);

   %%plot moving window:
   if data_type==2
      ylab1 = 'a, ms^{-2}';
      ylab2 = 'S, m^2s^{-3}';
   else
      ylab1 = 'w, m';
      ylab2 = 'S, m^2s';
   end


   for m=1:Nwindows
      t_start  = t_int(m,1);
      j0       = find(time==t_start);
      jj       = j0+(0:Nw-1)';
      time0    = time(jj);
      disp0    = displ(jj);
      disp0    = disp0-mean(disp0);
      %%
      figure(102);
      subplot(1,3,1),plot(time,displ,'k');
      Y  = 1.05*max(displ)*[-1 1];
      ylim(Y);
      hold on;
      plot(time0(1)+0*Y,Y,'r');
      plot(time0(end)+0*Y,Y,'r');
      hold off;
      GEN_proc_fig('time, s',ylab1)
      %%
      subplot(1,3,2);
      plot(time0,disp0,'k');
      ylim(Y);
      xlim([time0(1) time0(end)]);
      GEN_proc_fig('time, s',ylab1)
      ttl   = title(tstr);
      GEN_font(ttl);
      %%
      [Sn,ff]     = citeph_an2Sn(an_mat(:,m),df);
      [Smax,jmax] = max(Sn);
      periods     = 1./ff;
      Tp          = periods(jmax);
      fp          = ff(jmax);
      Sn_mat(:,m) = Sn;
      %periods(jmax-1:jmax+1)'

      %Sn_max   = max(max(Sn_mat));
      %SS       = Sn_max*1.10*[0 1];

      if DO_PLOT==2
         subplot(1,3,3);
         hold off;
         plot(T_target+0*SS,SS,'g')
         hold on;
         plot(periods,Sn,'k');
         xlim([0 4*T_target]);
         ylim(SS);
         GEN_proc_fig('period, s',ylab2)
         drawnow;
         hold off;
         pause(0.05);
      end
   end

   [PP,TM_vec] = meshgrid(1./ff,tm_vec);
   %{FF,TM_vec,Sn_mat}
   figure(111);
   %contour(TM_vec,PP,Sn_mat');
   pcolor(TM_vec,PP,Sn_mat');
   shading flat;
   ylim([0 5*T_target]);
   colorbar;
   ttl   = title(tstr);
   GEN_font(ttl);
   GEN_proc_fig('Time, s','Period, s');
   figname              = matfile;
   figname(end-2:end)   = 'eps';
   saveas(gcf,figname2,'epsc');
   disp(['saved: ' figname2]);
end
