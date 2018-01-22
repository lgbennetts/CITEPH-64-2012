function r_show_results(nbF)
% R_SHOW_RESULTS - Produce figures and movies showing the equilibrium
% position of floes over the entire domain, enabling validation of the
% data. This function can be expanded to do some other calculations and
% figures.
%
% Syntax:  r_show_results(nbF)
%
% Inputs:
%    nbF       - The number of floes (40 or 80)
%
% Example:
%    r_show_results(80)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: data_$nbF.mat (produced by R_COMBINE_RESULTS)
%
% See also: R_COMBINE_RESULTS, R_LIST, R_DATA_POS, R_DATA_POS_ALL
%
% Author: Dany Dumont
% UQAR/ISMER
% email: dany_dumont@uqar.ca
% Website: http://www.ismer.ca/dumont-dany
% November 2017
% ______________________________________________________________________

nbF = 80;
nbFStr  = num2str(nbF);
%% Open files
% rootdir   = r_root(nbF);        % You can change that
% datafile = [rootdir,'/data_',nbFStr,'f.mat'];
datafile = ['data_',nbFStr,'f.mat'];
load(datafile)
infofile = ['info.mat'];
load(infofile)

fps      = info.fps;
ave      = info.average;
fpsStr   = num2str(fps);
aveStr   = num2str(ave);

nTests       = length(S);
nRows        = 16;
nFloesPerRow = 5;

for i = 1:nTests
    
    hs      = S(i).hs./100;
    tp      = S(i).tp./10;
    expStr  = num2str(S(i).exp);
    nFrames = S(i).nFrames;
    x       = S(i).x;
    y       = S(i).y;
    xFloe   = S(i).xFloe;
    xRow    = S(i).xRow;
    yf_avg  = S(i).yf_avg;
    yi_avg  = S(i).yi_avg;
    if nbF == 80
        yRaft   = S(i).yRaft;
        rafted  = S(i).rafted;
        raftIntf = S(i).raftIntf;
        floeDist = S(i).floeDist;
    end
    
    % Here we can recompute the rafted array using another distance
    % threshold. A different value can be used to assess collisions instead
    % of rafting, based on a closer inspection of contact duration.
    cutoff = 0.96;
    raftIntf = floeDist;
    raftIntf(floeDist > cutoff)  = 0;
    raftIntf(floeDist <= cutoff) = 1;
    for n = 1:nFrames
        for r = 1:nRows
            tmp1 = squeeze(raftIntf(n,:,r));
            tmp2 = inv_and(tmp1);
%             rafted(n,:,r) = tmp2;
            rafted2(n,:,r) = tmp2;
        end
    end
    
    %% Make movie of individual tests showing rafting (and collisions?)
    
    if i == 10
        % Making a movie takes some time. It uses the VideoWriter function
        % which requires that you don't hide the figure frame during
        % recording. Uncomment the following relevant lines to do the
        % movie (see also at the end of the loop).
        v = VideoWriter([rootdir, ...
            '/miz_',fpsStr,'fps_',nbFStr,'f_',aveStr,'_',expStr,'.avi']);
        open(v);
        
        x = S(i).x;
        y = S(i).y;
        
        for n = 1:nFrames
            
            disp(['   ',num2str(n),' of ',num2str(nFrames)])
            figure(2); clf; hold on
            
            for r = 1:nRows
                for k = 1:nFloesPerRow
                    if n <= size(x,1)
                        if nbF == 80 && r == 7 || r == 10
                            circles(x(n,k,r),y(n,k,r),0.49,'EdgeColor','none', ...
                                'FaceColor',[0.8 0.2 0.2],'FaceAlpha',0.3)
                        elseif nbF == 40 && r == 6 || r == 10
                            circles(x(n,k,r),y(n,k,r),0.49,'EdgeColor','none', ...
                                'FaceColor',[0.8 0.2 0.2],'FaceAlpha',0.3)
                        else
                            if nbF == 80 && rafted(n,k,r)
                                circles(x(n,k,r),y(n,k,r),0.49,'EdgeColor','none', ...
                                    'FaceColor',[0.0 0.6 0.8],'FaceAlpha',0.3)
                            else
                                circles(x(n,k,r),y(n,k,r),0.49,'EdgeColor','none', ...
                                    'FaceColor',[0.0 0.1 0.8],'FaceAlpha',0.3)
                            end
                        end
                    end
                end
            end
            
            xlim([0 16])
            ylim([-6 2])
            text(1,1,['\itH_s=\rm',num2str(hs),' m \it T_p=\rm', ...
                num2str(tp),' s'], ...
                'FontName','Times New Roman','FontSize',16);
            text(13,-5.5,['\itt=\rm',num2str(n/5),' s'], ...
                'FontName','Times New Roman','FontSize',16);
            xlabel('\it x \rm (m)','FontName','Times New Roman','FontSize',18)
            ylabel('\it y \rm (m)','FontName','Times New Roman','FontSize',18)
            set(gca,'DataAspectRatio',[1 1 1],'Box','on', ...
                'FontName','Times New Roman','FontSize',16)
            hold off

            frame = getframe(gcf);
            writeVideo(v,frame);
        end
        close(v);
    end
    
    %% Makes figures of the equilibrim MIZ position
    
    figure(1); clf
    hold on
    g = 0.5.*[1 1 1];
    
    for r = 1:nRows
        for k = 1:nFloesPerRow
            if nbF == 80 && r == 7 || r == 10
                circles(xFloe(k,r),yf_avg(k,r),0.49,'EdgeColor','none', ...
                    'FaceColor',[0.8 0.2 0.2],'FaceAlpha',0.3)
                circles(xFloe(k,r),yi_avg(k,r),0.49, ...
                    'EdgeColor',g,'FaceColor','none')
            elseif nbF == 40 && r == 6 || r == 10
                circles(xFloe(k,r),yf_avg(k,r),0.49,'EdgeColor','none', ...
                    'FaceColor',[0.8 0.2 0.2],'FaceAlpha',0.3)
                circles(xFloe(k,r),yi_avg(k,r),0.49, ...
                    'EdgeColor',g,'FaceColor','none')
            else
                circles(xFloe(k,r),yf_avg(k,r),0.49,'EdgeColor','none', ...
                    'FaceColor',[0.0 0.1 0.8],'FaceAlpha',0.3)
                circles(xFloe(k,r),yi_avg(k,r),0.49, ...
                    'EdgeColor',g,'FaceColor','none')
            end
        end
    end
    plot([0 16],(mean(yf_avg(1,  2:2:16)) + 0.49).*[1 1],'--b')
    plot([0 16],(mean(yf_avg(end,2:2:16)) - 0.49).*[1 1],'--b')
    plot([0 16],(mean(yi_avg(1,  2:2:16)) + 0.49).*[1 1],'--','Color',g)
    plot([0 16],(mean(yi_avg(end,2:2:16)) - 0.49).*[1 1],'--','Color',g)
    
    xlim([0 16])
    ylim([-6 2])
    text(1,1,['\itH_s=\rm',num2str(hs),' m \it T_p=\rm', ...
        num2str(tp),' s'], ...
        'FontName','Times New Roman','FontSize',16);
    xlabel('\it x \rm (m)','FontName','Times New Roman','FontSize',18)
    ylabel('\it y \rm (m)','FontName','Times New Roman','FontSize',18)
    set(gca,'DataAspectRatio',[1 1 1],'Box','on', ...
        'FontName','Times New Roman','FontSize',16)
    hold off
    
    %pause(1)
    print([rootdir,'/mizeq_',fpsStr,'fps_',nbFStr,'f_',aveStr,'_',expStr,'.png'],'-dpng') ;
    
end

