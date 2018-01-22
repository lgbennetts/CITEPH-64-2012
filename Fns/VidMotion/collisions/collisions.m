function collisions

firstFrame = [ 607 2857 4683 5907 8157  9882 11757 13757 15582 18082];
lastFrame  = [1182 3182 5107 6732 8957 10682 12332 14232 16307 18507];

startTime  = firstFrame./25;
endTime    = lastFrame./25;

%% Undistortion step (comment if already done)

% for n=1:length(firstFrame)
%     r_undistort(201307311605,1,firstFrame(n),1,lastFrame(n));
% end

%% Interpolation step (comment if already done)

% for n=1:length(firstFrame)
%     %done! r_interp(201307311605,1,1,5,0);
%     % r_interp(201307311605,1,2,5,0);
%     r_interp(201307311605,1,2315+3,5,0);
%     % r_interp(201307311605,1,4,5,0);
% end
% 
% return

%% Collision parameters

nObjects = 2;
nCollisions = length(firstFrame);

load data_201307311605_1

% Sort floes by color

nFrames = length(nfloes);

tmp1 = nan.*ones(nFrames,nObjects);
tmp2 = nan.*ones(nFrames,nObjects);

for n = 1:nFrames

    for k = 1:nObjects
        if k == 1
            kind = find(color(n,:) == 3,1);
            if ~isempty(kind)
                tmp1(n,k) = xpos(n,kind);
                tmp2(n,k) = ypos(n,kind);
            end
        elseif k == 2
            kind = find(color(n,:) == 2,1);
            if ~isempty(kind)
                tmp1(n,k) = xpos(n,kind);
                tmp2(n,k) = ypos(n,kind);
            end
        else
        end
    end
end

xpos =  0.01.*tmp1 - 0.1;
ypos = -0.01.*tmp2 + 0.1;

for c = 1:nCollisions
    
    cstr = num2str(c);
    
    eval(['x_',cstr,'=xpos(time >= startTime(c) & time <= endTime(c),:);']);
    eval(['y_',cstr,'=ypos(time >= startTime(c) & time <= endTime(c),:);']);
    eval(['t1_',cstr,'=time(time >= startTime(c) & time <= endTime(c),:);']);
    
    x  = xpos(time >= startTime(c) & time <= endTime(c),:);
    y  = ypos(time >= startTime(c) & time <= endTime(c),:);
    t1 = time(time >= startTime(c) & time <= endTime(c),:);
    nFrames = length(t1);
    
    u  = nan.*ones(nFrames-1,nObjects);
    v  = nan.*ones(nFrames-1,nObjects);
    t2 = nan.*ones(nFrames-1,1);
    
    for n = 1:nFrames-1
        
        dt    = t1(n+1) - t1(n);
        t2(n) = t1(n) + 0.5.*dt;
        for k = 1:nObjects
            u(n,k) = (x(n+1,k) - x(n,k))./dt;
            v(n,k) = (y(n+1,k) - y(n,k))./dt;
        end
    end
    
    d     = sqrt((x(:,1) - x(:,2)).^2 + (y(:,1) - y(:,2)).^2);
    ic    = find(d == min(d),1,'first');
    dmin(c) = d(ic);
    tc(c) = t1(ic);
    
    dind = 22;
    
    x  = squeeze(x(ic-dind:ic+dind,:));
    y  = squeeze(y(ic-dind:ic+dind,:));
    u  = squeeze(u(ic-dind:ic+dind,:));
    v  = squeeze(v(ic-dind:ic+dind,:));
    d  = squeeze(d(ic-dind:ic+dind,:));
    t1 = squeeze(t1(ic-dind:ic+dind,:));
    t2 = squeeze(t2(ic-dind:ic+dind,:));
    
    ic = find(d == min(d),1,'first');
    
    xx1 = detrend(x(:,1),'linear',ic);
    xx2 = detrend(x(:,2),'linear',ic);
    yy1 = detrend(y(:,1),'linear',ic);
    yy2 = detrend(y(:,2),'linear',ic);
    dd  = detrend(d,'linear',ic);
    
    xt1 = x(:,1) - xx1;
    xt2 = x(:,2) - xx2;
    yt1 = y(:,1) - yy1;
    yt2 = y(:,2) - yy2;
    
    
    u1i(c) = (xt1(ic-1) - xt1(1))   ./(t1(ic-1) - t1(1));
    u1f(c) = (xt1(end)  - xt1(ic+1))./(t1(end)  - t1(ic+1));
    v1i(c) = (yt1(ic-1) - yt1(1))   ./(t1(ic-1) - t1(1));
    v1f(c) = (yt1(end)  - yt1(ic+1))./(t1(end)  - t1(ic+1));
    
    u2i(c) = (xt2(ic-1) - xt2(1))   ./(t1(ic-1) - t1(1));
    u2f(c) = (xt2(end)  - xt2(ic+1))./(t1(end)  - t1(ic+1));
    v2i(c) = (yt2(ic-1) - yt2(1))   ./(t1(ic-1) - t1(1));
    v2f(c) = (yt2(end)  - yt2(ic+1))./(t1(end)  - t1(ic+1));
    
%% Figures
    
    lw = 1.5;

    % Speed of floes
    figure(1)
    
    subplot(2,1,1)
    plot(t2,u,'o'); hold on
    stem(tc(c), 10,'--k');
    stem(tc(c),-10,'--k');
    plot(t2(1:ic-1),u1i(c).*ones(size(t2(1:ic-1))),'-r','LineWidth',lw);
    plot(t2(1:ic-1),u2i(c).*ones(size(t2(1:ic-1))),'-r','LineWidth',lw);
    plot(t2(ic+1:end),u1f(c).*ones(size(t2(ic+1:end))),'-r','LineWidth',lw);
    plot(t2(ic+1:end),u2f(c).*ones(size(t2(ic+1:end))),'-r','LineWidth',lw);
    plot(t2,zeros(size(t2)),'--k'); hold off
    set(gca,'FontName','Times','FontSize',18);
    ylim([-0.2 0.2])
    xlim([tc(c)-1 tc(c)+1])
    title(['Collision number ',cstr])
    ylabel('\itu \rm(m/s)')
    
    subplot(2,1,2)
    plot(t2,v,'o'); hold on
    stem(tc(c), 10,'--k');
    stem(tc(c),-10,'--k');
    plot(t2(1:ic-1),v1i(c).*ones(size(t2(1:ic-1))),'-r','LineWidth',lw);
    plot(t2(1:ic-1),v2i(c).*ones(size(t2(1:ic-1))),'-r','LineWidth',lw);
    plot(t2(ic+1:end),v1f(c).*ones(size(t2(ic+1:end))),'-r','LineWidth',lw);
    plot(t2(ic+1:end),v2f(c).*ones(size(t2(ic+1:end))),'-r','LineWidth',lw);
    plot(t2,zeros(size(t2)),'--k'); hold off
    set(gca,'FontName','Times','FontSize',18);
    ylim([-0.4 0.1])
    xlim([tc(c)-1 tc(c)+1])
    ylabel('\itv \rm(m/s)')
    xlabel('\itt \rm(s)')
    
    print(['uv_',cstr],'-dpng');
    
    % Position of floes
    figure(2)
    
    subplot(2,1,1)
    plot(t1,x,'o'); hold on
    plot(t1,x(:,1)-xx1,'-r','LineWidth',lw);
    plot(t1,x(:,2)-xx2,'-r','LineWidth',lw);
%     stem(collision(c), 10,'--k');
%     stem(collision(c),-10,'--k');
    hold off
    set(gca,'FontName','Times','FontSize',18);
    %ylim([-0.2 0.2])
    xlim([tc(c)-1 tc(c)+1])
    title(['Collision number ',cstr])
    ylabel('\itx \rm(m)')
    
    subplot(2,1,2)
    plot(t1,y,'o'); hold on
    plot(t1,y(:,1)-yy1,'-r','LineWidth',lw);
    plot(t1,y(:,2)-yy2,'-r','LineWidth',lw);
    stem(tc(c), 10,'--k');
    stem(tc(c),-10,'--k');
    hold off
    set(gca,'FontName','Times','FontSize',18);
    ylim([-4 -1])
    xlim([tc(c)-1 tc(c)+1])
    ylabel('\ity \rm(m)')
    xlabel('\itt \rm(s)')
    
    print(['xy_',cstr],'-dpng');
    
    % Trajectories of floes
    figure(3)
    plot(x,y,'o');
    set(gca,'FontName','Times','FontSize',18);
    xlim([ 0.0 4.0])
    ylim([-5.0 1.0])
    title(['Collision number ',cstr])
    ylabel('\ity \rm(m)')
    xlabel('\itx \rm(m)')
    
    print(['traj_',cstr],'-dpng');
    
    % Distance between floes
    figure(4)
    plot(t1,d,'o'); hold on
    plot(t1,d(:,1)-dd,'-r','LineWidth',lw);
    stem(tc(c), 1.1,'--k','BaseValue',0.99);
    hold off
    set(gca,'FontName','Times','FontSize',18);
    text(tc(c)-0.9,1.01,['\itt_c \rm= ',num2str(tc(c)),' s'],'FontName','Times','FontSize',18)
    xlim([tc(c)-1 tc(c)+1])
    ylim([0.95 1.25])
    title(['Collision number ',cstr])
    ylabel('\itd \rm(m)')
    xlabel('\itt \rm(s)')
    
    print(['dist_',cstr],'-dpng');


end

save('collisions','dmin','u1i','u1f','v1i','v1f','u2i','u2f','v2i','v2f');

