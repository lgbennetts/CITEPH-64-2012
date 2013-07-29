%% citeph_plot_movingFFT.m
%% Author: Timothy Williams
%% Date:   20130729, 10:43:04 CEST
function citeph_plot_movingFFT(test_type,test_num,sens_name)

dir0  = '/work/timill/CITEPH-data/results_preliminary/conc_79/regular_processed/';
if strcmp(test_type,'calib')
   [dirname,T_target,H_target,type,expt_name] = citeph_get_calib_prams(test_num);%%NB full scale
elseif strcmp(test_type,'c79')
   [dirname,T_target,H_target,type,expt_name] = citeph_get_c79_prams(test_num)%%NB full scale
else%%'c39'
   [dirname,T_target,H_target,type,expt_name] = citeph_get_c39_prams(test_num);%%NB full scale
end

if strcmp(sens_name(1),'A')
   dir1     = [dir0 '/' dirname '/' sens_name([1 3])];
elseif strcmp(sens_name(1),'S')
   dir1     = [dir0 '/' dirname '/' 'S'];
elseif strcmp(sens_name(end-3:end),'zoom')
   dir1     = [dir0 '/' dirname '/' 'S_zoom'];
end
matfile  = [dir1 '/' 'movingFFT_' sens_name];

load(matfile);

Nwindows = size(an_mat,2);

if DO_PLOT
   %%plot time series:
   figure(101);
   subplot(2,1,1);
   plot(tm_vec,Tp_vec,'k');
   hold on;
   plot(tm_vec,Tm02_vec,'--r');
   hold off;
   GEN_proc_fig('time, s','T_p & T_m02, s');
   tstr  = [expt_name,', ' sensor_name];
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

   %%plot moving window:
   if data_type==2
      ylab1 = 'a, ms^{-2}';
      ylab2 = 'S, m^2s^{-3}';
   else
      ylab1 = 'w, m';
      ylab2 = 'S, m^2s';
   end


   for m=m0:Nwindows
      jj          = (1+(m-1)*n_shift)+(0:Nw-1)';
      time0       = time(jj);
      disp0       = displ(jj);
      disp0       = disp0-mean(disp0);
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
      subplot(1,3,3);
      [Smax,jmax] = max(Sn_mat(:,m));
      periods  = 1./ff;
      Tp = periods(jmax);
      fp = ff(jmax);
      %periods(jmax-1:jmax+1)'

      Sn_max   = max(max(Sn_mat));
      SS       = Sn_max*1.10*[0 1];

      hold off;
      plot(T_target+0*SS,SS,'g')
      hold on;
      plot(periods,Sn_mat(:,m),'k');
      xlim([0 4*T_target]);
      ylim(SS);
      GEN_proc_fig('period, s',ylab2)
      drawnow;
      hold off;
      pause(0.05);
   end
end


function [ff,Sn] = an2Sn(an,df);

N  = length(an);
ff = (0:N/2-1)*df;
Sn = 0*ff;

Sn(1)    = abs(an(1))^2/df;
jj       = 2:N/2;
Sn(jj)   = abs( an(jj).^2 )/df + abs( an(N+1-jj).^2 )/df;
