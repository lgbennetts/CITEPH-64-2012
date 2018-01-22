function analysis
% This code performs the following:
%  1. plots drift of MIZ as a function of wave steepness
%  2. plots spread of MIZ as a function of wavelength and wave steepness
%  3. plots no. of events (contact, rafting & non-rafting collisions) as a
%      function of wave properties
%
% Lucas Yiew
% December 2017

close all

% addpath('LucasYiew')

load('data_80f.mat')

for i = 1:length(S)
 hs = S(i).hs./100; % [m]
 tp = S(i).tp./10;  % [s]
 
 [field] = wavefield('T',tp,3.1);
 k = field{5,2};
 ka(i) = k*hs/2;
 lambda(i) = 2*pi/k;
 T(i) = tp;
 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % CALCULATE DRIFT USING MEAN SPREAD
 for ii = 1:length(S(i).dy_avg(:,1))
  dy_avg(ii) = mean((S(i).dy_avg(ii,[1:3,5:14,16]))); % mean spread for each row
 end
 Dy_s_miz(i) = mean(dy_avg); % drift: mean spread for miz
%  Dy_s_miz(i) = mean(dy_avg)./(2*pi/k); % drift: mean spread for miz normalised wrt wavelength
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % CALCULATE DRIFT USING CENTRE DISPLACEMENT
 y_i = (mean(S(i).yi_avg(1,[1:3,5:14,16])) + mean(S(i).yi_avg(end,[1:3,5:14,16])))/2; % mean initial y of miz centre
 y_f = (mean(S(i).yf_avg(1,[1:3,5:14,16])) + mean(S(i).yf_avg(end,[1:3,5:14,16])))/2; % mean final y of miz centre
 Dy_d_miz(i) = y_f - y_i; % drift: mean centre of miz
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % CALCULATE OVERALL SPREAD OF MIZ
 Dy_i = (mean(S(i).yi_avg(1,[1:3,5:14,16])) - mean(S(i).yi_avg(end,[1:3,5:14,16]))) + 2*0.49; % mean initial overall spread
 Dy_f = (mean(S(i).yf_avg(1,[1:3,5:14,16])) - mean(S(i).yf_avg(end,[1:3,5:14,16]))) + 2*0.49; % mean final overall spread
%  dy_miz(i) = (Dy_f - Dy_i); % overall spread for miz 
%  dy_miz(i) = (Dy_f - Dy_i)./(2*pi/k); % overall spread for miz wrt wavelength
 dy_miz(i) = (Dy_f - Dy_i)./(hs/2); % overall spread for miz wrt amplitude
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % CALCULATE NUMBER OF COLLISIONS
 contact = S(i).rafted; % THIS IS DIFFERENT TO rafted BELOW ###
 % DANY'S CODE
%  yRaft   = S(i).yRaft;
%  rafted  = S(i).rafted;
 raftIntf = S(i).raftIntf;
 floeDist = S(i).floeDist;
 % Here we can recompute the rafted array using another distance
 % threshold. A different value can be used to assess collisions instead
 % of rafting, based on a closer inspection of contact duration.
 cutoff = 0.96;
 raftIntf = floeDist;
 raftIntf(floeDist > cutoff)  = 0;
 raftIntf(floeDist <= cutoff) = 1;
 nFrames = length(contact);
 for n = 1:nFrames
     for r = 1:16
         tmp1 = squeeze(raftIntf(n,:,r));
         tmp2 = inv_and(tmp1);
         rafted(n,:,r) = tmp2; % ###
     end
 end
 contact = rafted;
 % change [nFrames,5,16] array into [nFrames,4,16] array, showing only 
 % contact flags between each row.
 for j = 1:16
  for jj = 1:4
   for jjj = 1:length(contact)
    if contact(jjj,jj,j) == 1 && contact(jjj,jj+1,j) == 1
     dum(jjj,jj,j) = 1;
    elseif contact(jjj,jj,j) == 1 && contact(jjj,jj+1,j) == 2
     dum(jjj,jj,j) = 1;
    elseif contact(jjj,jj,j) == 2 && contact(jjj,jj+1,j) == 1
     dum(jjj,jj,j) = 1;
    elseif contact(jjj,jj,j) == 2 && contact(jjj,jj+1,j) == 2
     dum(jjj,jj,j) = 1;
    else
     dum(jjj,jj,j) = 0;
    end
   end
  end
 end
 clear contact; contact = dum; clear dum
 % calculate number of contact, rafting & nonrafting events
 nContacts = zeros(4,16);
 nRafting  = zeros(4,16);
 for j = 1:16
  for jj = 1:4
   if isnan(contact(:,jj,j))
   else
    tContact = 0;
    for jjj = 1:length(contact)-1
     % identify number of contact events
     if contact(jjj,jj,j) == 1 && contact(jjj+1,jj,j) == 0 
      nContacts(jj,j) = nContacts(jj,j) + 1;
      tContact = 0;
     % time of each contact
     elseif contact(jjj,jj,j) == 1 && contact(jjj+1,jj,j) == 1
      tContact = tContact + 1;
      % rafting duration defined here
      if tContact == 5 % sampling frequency = 5Hz therefore tContact = 5 => 1sec
       nRafting(jj,j) = nRafting(jj,j) + 1;
      end
     end
    end
   end
  end
 end
 % disregard columns 7 and 10
 for jj = 1:4
  nContacts(jj,7)  = 0;
  nContacts(jj,10) = 0;
  nRafting(jj,7)   = 0;
  nRafting(jj,10)  = 0;
 end
 clear nNRafting
 nNRafting = nContacts - nRafting;
 nEvents(i).Contacts = nContacts;
 nEvents(i).Rafting  = nRafting;
 nEvents(i).NRafting = nNRafting;
 nEvents(i).TContacts = sum(sum(nContacts));
 nEvents(i).TRafting  = sum(sum(nRafting));
 nEvents(i).TNRafting = sum(sum(nNRafting));
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT FIGURES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% MIZ DRIFT VS WAVE STEEPNESS
f1 = figure(1);
hold on
box on
plot(ka,Dy_s_miz,'kx')
plot(ka,Dy_d_miz,'b^')
xlabel('Wave Steepness, ka')
ylabel('MIZ Drift, \Deltay_{MIZ} [m]')
legend('Using Mean (Relative) Spread wrt Each MIZ Row','Using Mean MIZ-Centre Displacement')

% MIZ SPREAD VS WAVELENGTH & WAVE STEEPNESS
f2 = figure(2);
fs1= subplot(2,1,1);
box on
plot(lambda/(2*0.49),dy_miz,'rs')
xlabel('Wavelength / Floe Length')
% ylabel('MIZ Spread / Wave Amplitude, \deltay_{MIZ} / a ')
legend('Mean Overall Spread of MIZ')
%
fs2 = subplot(2,1,2);
box on
plot(ka,dy_miz,'rs')
xlabel('Wave Steepness, ka')
% ylabel('MIZ Spread / Wave Amplitude, \deltay_{MIZ} / a ')
legend('Mean Overall Spread of MIZ')
p1=get(fs1,'position');
p2=get(fs2,'position');
height=p1(2)+p1(4)-p2(2);
h3=axes('position',[p2(1) p2(2) p2(3) height],'visible','off');
h_label=ylabel('MIZ Spread / Wave Amplitude, \deltay_{MIZ} / a ','visible','on');

% NUMBER OF COLLISIONS AND RAFTING EVENTS VS WAVE PERIOD & STEEPNESS
% large amplitude tests S(3), S(12), S(13)
for j = 1:length(nEvents)
 TContacts(j) = nEvents(j).TContacts;
 TRafting(j)  = nEvents(j).TRafting;
 TNRafting(j) = nEvents(j).TNRafting;
 if TContacts(j) == 0
  TContacts(j) = -100;
 end
 if TRafting(j) == 0
  TRafting(j) = -100;
 end
 if TNRafting(j) == 0
  TNRafting(j) = -100;
 end
end
%
f3 = figure(3);
set(gcf, 'Position', [500, 10, 600, 1400])
subplot(3,1,1)
hold on
box on
img = imread('fig_BenWil15-Transmission.png');
image('CData',img,'XData',[0.505 2.02],'YData',[2000 0])
plot(T,TContacts,'mo')
plot(T,TRafting,'r^')
plot(T,TNRafting,'bv')
plot(1,100,'wx')
for j = [3,12,13]
 plot(T(j),TContacts(j),'mo','LineWidth',1.5)
 plot(T(j),TRafting(j),'r^','LineWidth',1.5)
 plot(T(j),TNRafting(j),'bv','LineWidth',1.5)
end
xlabel('Wave Period [s]')
ylabel('Number of Events')
legend('Contact Events','Rafting Events','Non-Rafting Events',...
 'Large Amplitude Tests (Bold Symbols)',...
 'Location','East')
xlim([0.5 2.05])
ylim([-100 2000])
set(gca,'YTick',[-100,0:200:600],'YTickLabel',{'No Contact','0','200','400','600'})
box on
plot([0.5 2.05],[0 0],'k--')
%
subplot(3,1,3)
hold on
box on
plot(ka,TContacts,'mo')
plot(ka,TRafting,'r^')
plot(ka,TNRafting,'bv')
for j = [3,12,13]
 plot(ka(j),TContacts(j),'mo','LineWidth',1.5)
 plot(ka(j),TRafting(j),'r^','LineWidth',1.5)
 plot(ka(j),TNRafting(j),'bv','LineWidth',1.5)
end
xlabel('Wave Steepness, ka')
ylabel('Number of Events')
set(gca,'YTick',[-100,0:200:600],'YTickLabel',{'No Contact','0','200','400','600'})
xlim([0 0.14])
plot([0 0.14],[0 0],'k--')
%
subplot(3,1,2)
hold on
box on
plot(lambda/(0.49*2),TContacts,'mo')
plot(lambda/(0.49*2),TRafting,'r^')
plot(lambda/(0.49*2),TNRafting,'bv')
for j = [3,12,13]
 plot(lambda(j)/(0.49*2),TContacts(j),'mo','LineWidth',1.5)
 plot(lambda(j)/(0.49*2),TRafting(j),'r^','LineWidth',1.5)
 plot(lambda(j)/(0.49*2),TNRafting(j),'bv','LineWidth',1.5)
end
xlabel('Wavelength / Floe Length')
ylabel('Number of Events')
set(gca,'YTick',[-100,0:200:600],'YTickLabel',{'No Contact','0','200','400','600'})
xlim([0 7])
plot([0 7],[0 0],'k--')

pause

% save figures
print(f1,'fig_mizdrift_vs_ka.png','-dpng') ;
print(f2,'fig_mizspread_vs_wl_ka.png','-dpng') ;
print(f3,'fig_contactevents_vs_wl.png','-dpng') ;


% pause
