

%% IPSC, later IPSC, EPSC, later EPSC; cf. Excel sheet
indizes = [0	0	0	0
0	0	0	0
0	0	0	0
0	0	0	0
0	0	0	0
1	1	1	1
0	0	0	1
0	1	1	1
1	1	1	1
0	1	0	1
0	0	0	0
0	0	0	0
1	1	1	1
0	0	1	1
0	0	0	0
0	0	0	0
1	1	1	1
0	0	0	0
0	0	0	0
0	1	0	1
1	1	1	1
1	1	1	1
0	0	0	0
0	0	0	1
0	0	0	0
0	0	0	1
0	0	1	1
0	1	0	0
1	0	1	1
0	1	0	1
0	1	0	1
0	1	0	0
1	1	0	0
0	0	0	0
0	1	0	1
1	1	1	1
1	1	0	0
1	1	1	1
1	1	1	1
0	0	0	0
1	1	1	1
1	1	1	1
1	1	1	1
1	1	1	1
1	0	0	0
0	0	0	0
0	0	1	1
1	1	0	1
0	0	0	0
1	1	1	1
0	0	0	0
0	0	0	0
0	0	1	1
0	0	0	1
0	0	1	1
1	1	0	1
1	1	1	1
0	1	0	1
0	0	0	0
0	0	0	0
1	1	0	1
1	1	1	1
1	1	1	1
1	1	1	1
1	1	1	1
1	1	1	1
1	1	1	1
1	1	1	1
1	1	1	1
1	1	1	1
1	1	0	1
1	1	1	1
1	1	1	1
1	1	1	1
1	1	1	1
1	1	1	1
1	1	1	1
1	1	1	1
1	1	1	1
1	1	1	1
1	1	1	1
1	1	1	1
1	1	1	1
1	1	1	1
0	1	0	1
1	1	1	1
1	1	1	1
1	1	0	0
1	1	1	1
1	1	1	1
0	0	1	1
1	1	1	1
1	1	1	1
1	1	1	1
0	0	0	0
1	1	0	1
1	1	1	1
1	1	1	1
0	0	1	1
1	1	1	1
1	1	1	1
1	1	0	1
1	1	0	0
0	0	0	0
0	0	0	0
0	0	0	1
1	1	1	1
1	1	1	1
1	1	1	1
1	1	1	1
1	1	1	1
1	1	0	0
0	0	1	1
1	1	1	1
1	1	1	1
1	1	1	1
1	1	1	1
1	1	1	1
0	1	0	1
0	1	0	1
0	1	0	1
0	1	0	1
1	1	1	1
1	1	0	1
1	1	0	1
1	1	0	1
1	1	0	1
0	1	0	1
0	1	0	1
0	1	0	1
0	1	0	1
1	1	0	1
1	1	1	1
0	1	0	1
1	1	1	1
1	1	1	1
1	1	1	1
1	1	1	1
1	1	0	1
1	1	0	1
1	1	0	0
1	1	0	0
0	0	0	1
0	0	1	1
0	1	0	1
0	0	1	1
0	0	0	1];

%% create onset dataset
YX = find(sum(indizes'));
kkkk = kkkk + 1; 
indizes(YX(kkkk),:)
YX(kkkk)
for kkkk = 77:123
    kkkk
    [onsetVC70{kkkk}, onsetVC00{kkkk}] = onsetSelect(YX(kkkk),[]);
end
save('tempOnset.mat','onsetVC70','onsetVC00');

%% create dataset of delays between EPSCs and IPSCs (AA) and variability (sAA) from manual annotation
counter = 1; nennersum = 0; totalsum = 0;
cmap = distinguishable_colors(3);
for j = 1:numel(onsetVC70)
    revIndex(counter) = YX(j);
    dateXnew(j) = str2double(datasetSingleCells{revIndex(counter)}.CellID(1:6));

end
dateUnique = unique(dateXnew); dateUnique = sort(dateUnique);

for j = 1:numel(onsetVC70)
    A = mean(onsetVC70{j},2);
    B = mean(onsetVC00{j},2);
    A(A<8) = NaN;
    B(B<8) = NaN;
    sA = std(onsetVC70{j},0,2);
    sB = std(onsetVC00{j},0,2);
    for k = 1:numel(A)
        AA(counter) = A(k)-B(k);
        sAA(counter) = sqrt(sA(k)^2+sB(k)^2);
        revIndex(counter) = YX(j);
        dateXnew(j) = str2double(datasetSingleCells{revIndex(counter)}.CellID(1:6));
        dateIndex = strfind(dateUnique,dateXnew(j));
%         errorbar(j,AA(counter),sAA(counter));
%         if dateXnew
        counter = counter + 1;
        figure(869); hold on; plot(dateIndex,A(k),'o','MarkerSize',4/sqrt(0.1+sA(k)*10),'Color',cmap(k,:)); plot(dateIndex,B(k),'.','MarkerSize',10/sqrt(0.1+sB(k)*10),'Color',cmap(k,:))
%         figure(889); hold on; plot(counter,A(k),'kx','MarkerSize',9); plot(counter,B(k),'b.','MarkerSize',16)
    end
end
ax = gca;
set(ax,'XTick',1:numel(dateUnique))
clear dd
for k = 1:numel(dateUnique); dd{k} = mat2str(dateUnique(k)); end
set(ax,'XTickLabel',dd)
totalsum/nennersum


%% plot distribution of delay between EPSCs and IPSCs, sorted by delay / uncertainty
BB = AA(~isnan(AA));
sBB = sAA(~isnan(AA));
for j = 1:2
    if j == 1
        [~,XI] = sort((BB));
    else
        [~,XI] = sort(abs(sBB));
    end
    CC = BB(XI);
    sCC = sBB(XI);
    figure(59), subplot(1,2,j); errorbar(CC*1000,sCC*1000,'.k')
    hold on; plot(-10:150,(-10:150)*0,'k'); hold off;
    ylabel('delay EPSC vs. IPSC [ms]');
    xlabel('double-responding neurons, sorted')
end

%% show delay in spatial map
BB = AA(~isnan(AA) & sAA<0.5);
revIndexB = revIndex(~isnan(AA) & sAA<0.5);
sBB = sAA(~isnan(AA) & sAA<0.5);
t1 = 0.2; t2 = -0.2;
BB(BB>t1) = t1;
BB(BB<t2) = t2;
NeuronViewerMap(BB,revIndexB,72)






%% absolute odor onset for different fish and different odors
% manual selection of bad trials
for IX = 1:98
    clear onsetVC00 onsetVC70
    showRawDataOnlyVC(IX);
    num_dice = input('Enter trials to be discarded: ');
    [onsetVC70, onsetVC00,odorsVC,dateID] = extractVC(IX,num_dice);
    X_all{IX}.onsetVC70 = onsetVC70;
    X_all{IX}.onsetVC00 = onsetVC00;
    X_all{IX}.odorsVC = odorsVC;
    X_all{IX}.dateID = dateID;
end
save('GoodVCtrials_onsetAnalysis.mat','X_all','-v7.3')
% grouping according to fish
for IX = 1:98
    day{IX} = X_all{IX}.dateID;
end
days = unique(day);
for j = 1:numel(days)
    cells = find(~cellfun(@isempty,strfind(day,days{j})));
    for p = 1:3; vc70{p} = []; vc00{p} = []; end
    for ii = 1:numel(cells)
        for k = 1:numel(X_all{cells(1)}.odorsVC) % number of odors
            XX = squeeze(X_all{cells(ii)}.onsetVC70(k,:,:));
            if size(XX,1) < 1e4; XX = XX'; end
            vc70{k} = cat(2,vc70{k},XX);
            try;
                YY = squeeze(X_all{cells(ii)}.onsetVC00(k,:,:));
                if size(YY,1) < 1e4; YY = YY'; end
                vc00{k} = cat(2,vc00{k},YY);
            end
        end
    end
    figure(41), cmap = distinguishable_colors(3);
    for k = 1:3
        if ~isempty(vc00{k})
            Y1 = conv2(vc00{k},fspecial('gaussian',[10 1],7),'same'); Y1 = Y1';
            for i = 1:size(Y1,1); Y1(i,:) = (Y1(i,:) - mean(Y1(i,:)))/std(Y1(i,:)); end
            Y2 = conv2(vc70{k},fspecial('gaussian',[10 1],7),'same'); Y2 = Y2';
            for i = 1:size(Y2,1); Y2(i,:) = (Y2(i,:) - mean(Y2(i,:)))/std(Y2(i,:)); end
            subplot(2,1,1);
            plot(cumsum(nanmean(Y1)-mean(nanmean(Y1(:,1:0.8e5)))),'Color',cmap(1,:)); hold on;
            plot(cumsum(nanmean(Y2)-mean(nanmean(Y2(:,1:0.8e5)))),'Color',cmap(2,:)); hold off;
            xlim([7.5e4 15e4])
            subplot(2,1,2);
            plot(nanmean(Y1),'Color',cmap(1,:)); hold on;
            plot(nanmean(Y2),'Color',cmap(2,:)); hold off;
            xlim([7.5e4 15e4])
            ylim([-2 2]);
            xlabel(X_all{cells(1)}.dateID)
%             temp = ginput(1);
%             onset(j,k) = temp(1);
        end
    end
end
ordering = [1 NaN NaN;    1 NaN NaN;    1 NaN NaN;    1 NaN NaN;
    1 2 NaN;    1 2 NaN;    3 2 1;    3 2 1;
    3 2 1;    3 2 1;    3 2 1;    3 2 1;
    3 2 1;    3 2 1;    3 2 1;    1 2 3;
    1 2 3;    1 2 3;    2 1 3;    2 1 3;    2 1 3];
ordering(isnan(ordering)) = 3;
for p = 1:size(onset,1); onsetX(p,:) = onset(p,ordering(p,:)); end

onsetX(onsetX<8e4) = NaN;
for p = 1:size(onsetX,1)
    onsetX2(p,:) = onsetX(p,:) - nanmean(onsetX(p,:));
end
onset1 = nanmean(onsetX2(:,1));
onset2 = nanmean(onsetX2(:,2));
onset3 = nanmean(onsetX2(:,3));

for p = 1:size(onsetX,1)
    onsetZ(p,:) = nanmean(onsetX(p,:) - [onset1 onset2 onset3]) + [onset1 onset2 onset3];
end
onsetZ(5,:) = NaN;
for j = 1:3
    good_indizes = find(~isnan(onsetZ(:,j)));
    onsetZ(:,j) = interp1(good_indizes,onsetZ(good_indizes,j),1:size(onsetZ,1),'Nearest','extrap');
end

figure(68); plot(onsetZ), hold on; plot(onset,'.','Markersize',16)

save('OnsetDays.mat','onset','ordering','onsetZ')


%% 2D plot analog to Ming's timing plot; normalize with respect to noise before response onset
counter = 1;
FULL70 = zeros(300,2e5);
FULL00 = zeros(300,2e5);
for IX = 1:98
    IX
    onsetVC70 = X_all{IX}.onsetVC70; onsetVC70(onsetVC70 == 0) = NaN;
    onsetVC00 = X_all{IX}.onsetVC00; onsetVC00(onsetVC00 == 0) = NaN;
    odorsVC = X_all{IX}.odorsVC;
    for p = 1:numel(odorsVC)
        
        odorIndex = find(~cellfun(@isempty,strfind(datasetSingleCells{IX}.odors,odorsVC{p})));
        dayIndex = find(~cellfun(@isempty,strfind(days,X_all{IX}.dateID)));
        delay = onsetZ(dayIndex,odorIndex);

        onsetVC70(p,:,:) = circshift(onsetVC70(p,:,:),[0 0 1e5-round(delay)]);
        FULL70(counter,:) = (nanmean(squeeze(onsetVC70(p,:,:))) - nanmean(nanmean(squeeze(onsetVC70(p,:,1:0.8e5)))))/nanstd(nanmean(squeeze(onsetVC70(p,:,1:0.8e5))));
        try
            onsetVC00(p,:,:) = circshift(onsetVC00(p,:,:),[0 0 1e5-round(delay)]);
            FULL00(counter,:) = (nanmean(squeeze(onsetVC00(p,:,:))) - nanmean(nanmean(squeeze(onsetVC00(p,:,1:0.8e5)))))/nanstd(nanmean(squeeze(onsetVC00(p,:,1:0.8e5))));
        end
        counter = counter + 1;
    end
end
FULL70(FULL70 ==0) = NaN;
FULL00(FULL00 ==0) = NaN;
Crit1 = nanmean(FULL70(:,1e5:1.1e5),2);
Crit2 = nanmean(FULL00(:,1e5:1.05e5),2);
[ix, xi] = sort(Crit1);
% [ix, xi] = sort(Crit2);
Full70 = FULL70(xi,:);
Full00 = FULL00(xi,:);

figure(41); subplot(1,2,1); imagesc([0:20],1:300,Full70(:,1:10:end),[-15 3])
subplot(1,2,2);  imagesc([0:20],1:300,Full00(:,1:10:end),[0 30])

% difference (EPSC vs. IPSC) map in timing
DIFFull = Full70(:,1:10:end)+Full00(:,1:10:end);
Crit1 = nanmean(DIFFull(:,1e4:1.05e4),2);
% Crit2 = nanmean(DIFFull(:,1.025e4:1.05e4),2);
[ix, xi] = sort(Crit1);
DIFFullX = DIFFull(xi,:);
DIFFullX(~isfinite(DIFFullX)) = NaN;
figure(182), subplot(2,1,2); plot((1:20000)/1e3,smooth(nanmean(DIFFullX),10))
subplot(2,1,1); imagesc([0:20],1:300,DIFFullX,[-15 3])


%% quantify area-under-curve after supposed response onset in a 2 sec-time window
counter = 1;
FULL70 = zeros(300,2e5);
FULL00 = zeros(300,2e5);
for IX = 1:98
    IX
    onsetVC70 = X_all{IX}.onsetVC70; onsetVC70(onsetVC70 == 0) = NaN;
    onsetVC00 = X_all{IX}.onsetVC00; onsetVC00(onsetVC00 == 0) = NaN;
    odorsVC = X_all{IX}.odorsVC;
    for p = 1:numel(odorsVC)
        
        odorIndex = find(~cellfun(@isempty,strfind(datasetSingleCells{IX}.odors,odorsVC{p})));
        dayIndex = find(~cellfun(@isempty,strfind(days,X_all{IX}.dateID)));
        delay = onsetZ(dayIndex,odorIndex);

        onsetVC70(p,:,:) = circshift(onsetVC70(p,:,:),[0 0 1e5-round(delay)]);
        try
            onsetVC00(p,:,:) = circshift(onsetVC00(p,:,:),[0 0 1e5-round(delay)]);
        end
        for tt = 30
            t_end = 1e5+tt*500;
            t_start = 1e5+0*tt*500-500; % IMPORTANT : either integrate or moving window
            AUC70(p,IX,tt) = nanmean(nanmean(squeeze(onsetVC70(p,:,t_start:t_end))));
            AUC70std(p,IX,tt) = nanstd(squeeze(nanmean(onsetVC70(p,:,t_start:t_end),3))-squeeze(nanmean(onsetVC70(p,:,1:0.9e5),3)));
            try
                AUC00(p,IX,tt) = nanmean(nanmean(squeeze(onsetVC00(p,:,t_start:t_end))));
                AUC00std(p,IX,tt) = nanstd(squeeze(nanmean(onsetVC00(p,:,t_start:t_end),3))-squeeze(nanmean(onsetVC00(p,:,1:0.9e5),3)));
            catch
                disp(['Problem with cell no.',32,num2str(IX)]);
            end
        end
        nAUC70std(p,IX) = size(onsetVC70,2);
        pAUC70(p,IX) = nanmean(squeeze(nanmean(onsetVC70(p,:,1:0.9e5),3)));
        pAUC70std(p,IX) = nanstd(squeeze(nanmean(onsetVC70(p,:,1:0.9e5),3)));
        try
            nAUC00std(p,IX) = size(onsetVC00,2);
            pAUC00(p,IX) = nanmean(nanmean(squeeze(onsetVC00(p,:,1:0.9e5))));
        catch
            disp(['Problem with cell no.',32,num2str(IX)]);
        end
    end
end

%% plot average input integrators
timeX = (1:100)*50; % 100 refers to tt=1:100 as defined above
figure(47);
AUC70(AUC70==0) = NaN;
AUC00(AUC00==0) = NaN;
for rr = 1:size(AUC70,2)
    for pp = 1:size(AUC70,1)
        subplot(1,3,1); plot(timeX,-squeeze(AUC70(pp,rr,:)-pAUC70(pp,rr,:)),'k'); hold on;
        subplot(1,3,2); plot(timeX,squeeze(AUC00(pp,rr,:)-pAUC00(pp,rr,:)),'r'); hold on;
    end
end
normalize1 = sum(squeeze(nanmean(nanmean(AUC00-repmat(pAUC00,[1 1 100]),1),2)));
normalize2 = sum(squeeze(nanmean(nanmean(AUC70-repmat(pAUC70,[1 1 100]),1),2)));
subplot(1,3,3); plot(timeX,(squeeze(nanmean(nanmean(AUC70-repmat(pAUC70,[1 1 100]),1),2)))/normalize2,'k');
hold on; plot(timeX,squeeze(nanmean(nanmean(AUC00-repmat(pAUC00,[1 1 100]),1),2))/normalize1,'r');
xlabel('Time [msec]');
ylabel('Average current [pA]');
subplot(1,3,1); xlabel('Time [msec]'); ylabel('Average current [pA]'); hold off; ylim([-20 100])
subplot(1,3,2); xlabel('Time [msec]'); ylabel('Average current [pA]'); hold off;  ylim([-20 100])
subplot(1,3,3); hold off;

%% Timing comparison between IPSCs and EPSCs
counter = 0;
figure(49); 
for rr = 1:size(AUC70,2)
    if (indizes(rr,1) == 1 && indizes(rr,3) == 1)
        for pp = 1:size(AUC70,1)
            if sqrt(nAUC70std(pp,rr)*nAUC00std(pp,rr)) ~= 0 && ~isnan(sum(AUC00(pp,rr,:))) && ~isnan(sum(AUC70(pp,rr,:)))
                trajectoryx = squeeze(AUC70(pp,rr,:)); 
                trajectoryy = squeeze(AUC00(pp,rr,:));
                
                if 1% max(trajectoryy)-min(trajectoryy) > 30 || max(trajectoryx)-min(trajectoryx) > 30
                    counter = counter + 1;

                    trajectoryx = -trajectoryx; trajectoryx = trajectoryx - min(trajectoryx); exc = max(trajectoryx); trajectoryx = trajectoryx/max(trajectoryx);
                    trajectoryy = trajectoryy - min(trajectoryy); inh = max(trajectoryy); trajectoryy = trajectoryy/max(trajectoryy);

                    subplot(10,14,counter); 
                    plot(50,0,'MarkerSize',nthroot(inh,2)*7,'Color','r','Marker','.');
                    hold on; plot(100,0,'MarkerSize',nthroot(exc,2)*7,'Color','k','Marker','.');
                    plot(trajectoryx,'k'); hold on; plot(trajectoryy,'r'); hold off;
                    set(gca,'xtick',[],'ytick',[]);
                    set(gcf,'Color',[1 1 1]);
                    box off;

                    xlabel([datasetSingleCells{rr}.odors{pp},32,num2str(rr)],'Interpreter','None');
                end
            end
        end
    end
end

%% Timing comparison between IPSCs and EPSCs, histogram, take only strongly responding neurons
nn1 = zeros(100);
nn2 = zeros(100);
counter = 0
figure(547); cmap = distinguishable_colors(size(AUC70,2));
for rr = 1:size(AUC70,2)
    if (indizes(rr,1) == 1 && indizes(rr,3) == 1)
        for pp = 1:size(AUC70,1)
            
                trajectoryx = smooth(squeeze(AUC70(pp,rr,:)),1); 
                trajectoryy = smooth(squeeze(AUC00(pp,rr,:)),1);
    %         figure, plot(-trajectoryx,trajectoryy,'Color',cmap(1,:)); hold on;

            if max(trajectoryy)-min(trajectoryy) > 35 || max(trajectoryx)-min(trajectoryx) > 35
                counter = counter + 1
                [~,xx] = max((trajectoryy-min(trajectoryy)).^2+(trajectoryx-min(trajectoryx)).^2);

                trajectoryy = (trajectoryy - min(trajectoryy))/(max(trajectoryy) - min(trajectoryy));
                trajectoryx = (-trajectoryx - min(-trajectoryx))/(max(-trajectoryx) - min(-trajectoryx));
                figure(547); subplot(1,2,1); cmap = cbrewer('seq', 'YlOrRd', 20);
%                 plot(trajectoryx(:),trajectoryy(:),'Color',[0.8 0.8 0.8]); hold on;
                for k=1:30
                    plot(trajectoryx(k:k+1),trajectoryy(k:k+1),'Color',cmap(min(20,k),:),'Marker','.'); hold on;
                    plot(trajectoryx(k:k+1),trajectoryy(k:k+1),'Color',cmap(min(20,k),:)); hold on;
%                     plot(trajectoryx(xx+1:end),trajectoryy(xx+1:end),'Color',cmap(2,:),'Marker','.'); hold on;
                end
                figure(548); subplot(1,2,1); cmap = distinguishable_colors(2);
                plot(trajectoryx(1:xx),trajectoryy(1:xx),'Color',cmap(1,:),'Marker','.'); hold on;
                plot(trajectoryx(1:xx),trajectoryy(1:xx),'Color',cmap(1,:)); hold on;
                plot(trajectoryx(xx+1:end),trajectoryy(xx+1:end),'Color',cmap(2,:),'Marker','.'); hold on;
                plot(trajectoryx(xx:end),trajectoryy(xx:end),'Color',cmap(2,:)); hold on;
                nn1 =  nn1 + hist3([trajectoryx(1:xx), trajectoryy(1:xx)],[100 100]);
                nn2 = nn2 + hist3([trajectoryx(xx+1:end), trajectoryy(xx+1:end)],[100 100]);
                pause(0.001);

            end
        end
    end
end
axis([0 1 0 1]); drawnow; hold off;
subplot(1,2,2)
imagesc([0:0.1:1],[0:0.1:1],conv2(nn2',fspecial('gaussian',[8 8],5),'same'),[0 2.5]);
axis xy
figure(88); aa = zeros(22); aa(2) = 20; cmap = cbrewer('seq', 'YlOrRd', 20); imagesc(aa); colormap(cmap)


%% movie of trajectory in EPSC/IPSC space
% figure(41); ploterr(-(AUC70(:)-pAUC70(:)),AUC00(:)-pAUC00(:),pAUC70std(:)./nAUC70std(:),pAUC00std(:)./nAUC00std(:),'k.','logxy')
G = figure(44);
for tt = 1:100
    AUC70x = AUC70(:,:,tt);
    AUC00x = AUC00(:,:,tt);
    plot(-(AUC70x(:)-pAUC70(:)),AUC00x(:)-pAUC00(:),'k.','MarkerSize',14);
    axis([-5 80 -10 110]);
    title(['time = ',num2str(tt*50),32,'msec']);
    xlabel('average excitatory current, averaged over 200 msec [pA]')
    ylabel('average inhibitory current, averaged over 200 msec [pA]')
%     h = getframe(G);
%     imwrite(h.cdata,['frame',sprintf('%.3d',tt),'.png']);
    pause(0.1)
end


%% AUC over a certain period in time
figure(899);
t_end = 30; % 30 = 1.5 sec
AUC70x = AUC70(:,:,t_end);
AUC00x = AUC00(:,:,t_end);
cmap = jet(size(AUC70x,2));
for j = 1:size(AUC70x,2)
    for k = 1:size(AUC70x,1)
        
        xi = find(revIndex == j);
%         if isempty(xi) || k > numel(xi)
%             MSize(k,j) = 1;
%         else
%             MSize(k,j) = 1/(abs(AA(xi(k)))+0.02); 
%         end
        if ~isnan(MSize(k,j))
            nn = sqrt(nAUC70std(k,j)*nAUC00std(k,j));
            nn2 = min(nn/7,1);
            plot(-(AUC70x(k,j)-pAUC70(k,j))/20*t_end,(AUC00x(k,j)-pAUC00(k,j))/20*t_end,'.','Color',[1-nn2 1-nn2 1-nn2],'MarkerSize',max(1,nn*3));
            hold on;
        end
    end
end
xlabel('charge from EPSCs over 1.5 sec [pC]')
ylabel('charge from IPSCs over 1.5 sec [pC]')
hold off;



%% visualize parameters on the 2D projection map
t_end = 1; % 30 = 1.5 sec
AUC70x = AUC70(:,:,t_end)-pAUC70;
AUC00x = AUC00(:,:,t_end)-pAUC00;
propertyX = AUC00x+AUC70x;
% propertyX = AUC00x;
% propertyX(propertyX>quantile(propertyX(:),0.98)) = quantile(propertyX(:),0.98);
neuroIndex = [1:size(propertyX,2);1:size(propertyX,2);1:size(propertyX,2)];
ix = find(propertyX~=0 & ~isnan(propertyX));
propertyX = propertyX(ix); neuroIndex = neuroIndex(ix);
NeuronViewerMap(propertyX,neuroIndex,78)