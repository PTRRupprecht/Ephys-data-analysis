%% go to main folder




%% show electrophysiological recordings
FileList = dir('*AAAE*.xsg');
temp = zeros(500,2001,numel(FileList));
temp2 = zeros(500,2001,numel(FileList));
temp3 = zeros(500,2001,numel(FileList));
for i = 1:numel(FileList)
    i
    load(FileList(i).name,'-mat');
    samplerate = header.ephys.ephys.sampleRate;

    A = data.ephys.trace_2;
    L_trace = A;
    xdata = [1:size(L_trace)]/samplerate; %in sec
    new_data = timeseries(L_trace,xdata);
    idealfilter_data = idealfilter(new_data,[20 40],'pass');
    B = data.ephys.trace_1;
    if mean(B(:)) < 0; B = -B; trial_type = 'exc'; else trial_type = 'inh';  end % inverts inhibitory currents
    P_trace = B;
    xdata = [1:size(P_trace)]/samplerate; %in sec
    new_data = timeseries(P_trace,xdata);
    idealfilter_dataP = idealfilter(new_data,[0 100],'pass');
%     idealfilter_dataP0 = idealfilter(new_data,[0 500],'pass');
%     figure, plot(idealfilter_dataP)
%     figure(3), plot(idealfilter_dataP0)
    
    clear offsetk
    for k = 1:500
        offset = 110000+k*30;
        offsetk(k) = offset/1e4;
        X1 = idealfilter_dataP.Data(offset:offset+1000);
        X2 = idealfilter_dataP.Data(offset:offset+1000);%idealfilter_dataP.Data(offset:offset+1000);
        temp(k,:,i) = xcorr(X1,X2);
        X1 = idealfilter_data.Data(offset:offset+1000);
        X2 = idealfilter_data.Data(offset:offset+1000);%idealfilter_dataP.Data(offset:offset+1000);
        temp2(k,:,i) = xcorr(X1,X2);
        X1 = idealfilter_data.Data(offset:offset+1000);
        X2 = idealfilter_dataP.Data(offset:offset+1000);%idealfilter_dataP.Data(offset:offset+1000);
        temp3(k,:,i) = xcorr(X1,X2);
    end
    figure(31);
    subplot(ceil(sqrt(numel(FileList))),ceil(sqrt(numel(FileList))),i);
    imagesc([-1000:1000]/10,offsetk,temp(:,:,i))
    caxis([-1 1]*1e4)%3e4)
    figure(32);
    subplot(ceil(sqrt(numel(FileList))),ceil(sqrt(numel(FileList))),i);
    imagesc([-1000:1000]/10,offsetk,temp2(:,:,i))
    caxis([-1 1]*1e4)%*3e4)
    figure(33);
    subplot(ceil(sqrt(numel(FileList))),ceil(sqrt(numel(FileList))),i);
    imagesc([-1000:1000]/10,offsetk,temp3(:,:,i))
    caxis([-1 1]*1e4)%*3e4)

end
%% 21Dec2016
% AAAC: [3 4 5 10 11 12 17 18 19 23 24] [6 7 8 9 13 14 15 20 21 22] [200:250]
% AAAD: [3 4 5 10 11 12 13 20 21 22 23] [6 7 8 9 14 15 16 17 18 19] [220:260]
% AAAB: [3 4] [5 6 7 8] [240:290]
% AAAE: [4 5 6 7 12 13 14 15 20 21 22] [8 9 10 11 16 17 18 19] [250:300]
% AAAF: [1 2 3 4 5 10 11 12 13 18 19 20 21] [6 7 8 9 14 15 16 17 22 23] [260:300]
%% 14Dec2016
% AAAF: [1 2 3 4] [5 6 7 8]
% AAAG: [3 4 9 10 11 12 20 21 22 23 24 25 26 27] [5 6 7 8 13 14 15 16 17 18 19 28 29 30 31]
% AAAH: [2 3 9 10 11 12 13] [4 5 6 7 8]
%% 10Jan2017
% AAAF: [3 4 5 9 10 11 12 13 14 22 24 25 26] [6 7 8 15 16 17 18 19 21] [ 250:300] [11.7 12.0] [11000];
% AAAG1: [3 4 5 10 11 12 13 14 15] [6 7 8 9] [188:250] [] [11000];
% AAAG2: [23 24 25 29 30 31 32 34] [20 21 22 26 27 28] [188:250] [] [11000];

clear DD
for j = 1:200
    j
    indizes_exc = [4 5 6 7 12 13 14 15 20 21 22] ;
    indizes_inh = [8 9 10 11 16 17 18 19];
    indizes_exc = indizes_exc(randperm(length(indizes_exc)));
    indizes_inh =  indizes_inh(randperm(length(indizes_inh)));
    %mock 1
    indizes_exc1 = sort(indizes_inh(1:floor(numel(indizes_inh)/2)));
    indizes_inh1 =  sort(indizes_inh(floor(numel(indizes_inh)/2)+1:end));
    %mock2
    indizes_inh2 = indizes_exc(1:floor(numel(indizes_exc)/2));
    indizes_exc2 =  indizes_exc(floor(numel(indizes_exc)/2)+1:end);

%     LFP_exc = mean(temp(:,:,indizes_exc),3);
%     LFP_inh =  mean(temp(:,:,indizes_inh),3);
%     VC_exc = mean(temp2(:,:,indizes_exc),3);
%     VC_inh =  mean(temp2(:,:,indizes_inh),3);
    XX_exc = mean(temp3(:,:,indizes_exc),3);
    XX_inh =  mean(temp3(:,:,indizes_inh),3);
%   figure, imagesc(XX_exc); % figure, imagesc(XX_inh)
    XX_exc1 = mean(temp3(:,:,indizes_exc1),3);
    XX_inh1 =  mean(temp3(:,:,indizes_inh1),3);
    XX_exc2 = mean(temp3(:,:,indizes_exc2),3);
    XX_inh2 =  mean(temp3(:,:,indizes_inh2),3);


%     [~,IX] = max(mean(XX_exc(225:275,:),1));
%     [~,IX2] = max(mean(XX_inh(225:275,:),1));
%     [~,IX] = max(mean(XX_exc(75:200,:),1));
%     [~,IX2] = max(mean(XX_inh(75:200,:),1));

    rangge = 250:300;

    [~,IX] = max(mean(XX_exc(rangge,:),1));
    [~,IX2] = max(mean(XX_inh(rangge,:),1));
    (IX-IX2)/10
    [~,IX] = max(mean(XX_exc1(rangge,:),1));
    [~,IX2] = max(mean(XX_inh1(rangge,:),1));
    DD(j) = (IX-IX2)/10;
    [~,IX] = max(mean(XX_exc2(rangge,:),1));
    [~,IX2] = max(mean(XX_inh2(rangge,:),1));
    DD(j+200) = (IX-IX2)/10;
end
DD2 = reshape(DD,[100 4]); DD2 = mean(DD2,2);

std(DD2)

figure(999), plot((-1000:1000)/10,mean(XX_exc(rangge,:),1)/max(mean(XX_exc(rangge,:),1))); hold on; plot((-1000:1000)/10,mean(XX_inh(rangge,:),1)/max(mean(XX_inh(rangge,:),1)),'r');



figure(22)
subplot(2,3,1); imagesc(LFP_exc)
subplot(2,3,2); imagesc(LFP_inh)
subplot(2,3,3); imagesc(VC_exc)
subplot(2,3,4); imagesc(VC_inh)
subplot(2,3,5); plot((-1000:1000)/10,mean(XX_exc(250:350,:),1))
subplot(2,3,6); plot((-1000:1000)/10,mean(XX_inh(250:350,:),1))



% 21Dec2016
% AAAB : 3.0 % 1.055
% AAAC : 2.6 % 0.528
% AAAD : 3.1 % 0.548
% AAAE : 2.6 % 0.392
% AAAF : 3.6 % 0.868
% 14Dec2016 - bad quality, discarded
% AAAG : 0.5 % 0.59 % 0.80 > 0.70
% AAAH : -1.1 % 0.48 % 0.88 > 0.68
% 10Jan20017
% AAAF : 3.6 % 0.552
% AAAG1 : 4.0 % 0.698
% AAAG2 : 2.1 % 0.846
AllData = [3.0 2.6 3.1 2.6 3.6 3.6 4.0 2.1;
    1.055 0.528 0.548 0.392 0.868 0.552 0.698 0.846];

AllData = 

Xd = AllData(1,:);
Wd = 1./AllData(2,:).^2;

errror = sqrt(sum(AllData(2,:).^2))/numel(AllData(2,:))

% (2.99 +- 0.25) ms

reorder = randperm(size(AllData,2));
AllData = AllData(:,reorder);

figure(44),
plot([0 size(AllData,2)+1],[1 1]*2.99,'Color',[0.5 0.5 0.5]); hold on;
plot([0 size(AllData,2)+1],[1 1]*(2.99-2*errror),'--','Color',[0.5 0.5 0.5]); hold on;
plot([0 size(AllData,2)+1],[1 1]*(2.99+2*errror),'--','Color',[0.5 0.5 0.5]); hold on;
errorbar(1:size(AllData,2),AllData(1,:),AllData(2,:),'.k','MarkerSize',18,'LineWidth',2)
xlim([0.5 size(AllData,2)+0.5]); ylim([-10 10]);
hold on; plot([0 size(AllData,2)+1],[0 0],'k');
xlabel('Cell Index'); ylabel('Relative delay exc vs. inh [ms]');
ylim([-9 9])
hold off

figure(23); plot(-2000:2000,xcorr(mean(XX_exc(225:300,:),1),mean(XX_inh(225:300,:),1)))

clear selectedCycles
%% show electrophysiological recordings
FileList = dir('*AAAB*.xsg');
cmap = lines(numel(FileList));
for i = 1:numel(FileList)
    load(FileList(i).name,'-mat');
    samplerate = header.ephys.ephys.sampleRate;

    A = data.ephys.trace_2;
    L_trace = A;
    xdata = [1:size(L_trace)]/samplerate; %in sec
    new_data = timeseries(L_trace,xdata);
    idealfilter_data = idealfilter(new_data,[20 40],'pass');
%     idealfilter_data2 = idealfilter(new_data,[2 60],'pass');
%     idealfilter_data3 = idealfilter(new_data,[2 100],'pass');
%     idealfilter_data4 = idealfilter(new_data,[2 200],'pass');
    
    B = data.ephys.trace_1;
    if mean(B(:)) < 0; B = -B; trial_type = 'exc'; else trial_type = 'inh';  end % inverts inhibitory currents
    P_trace = B;
    xdata = [1:size(P_trace)]/samplerate; %in sec
    new_data = timeseries(P_trace,xdata);
    idealfilter_dataP = idealfilter(new_data,[20 40],'pass');
    
    global pos_manually
    pos_manually = [];
    figure(444), plot((1:numel(B))/1e4,B,'k'); xlim([11.3 13.5])
    figure(i), hold on;
    plot(idealfilter_data(1:10:end)-mean(idealfilter_data(1:10:end)),'Color',cmap(1,:),'LineWidth',2)
    plot(idealfilter_dataP(1:10:end)-mean(idealfilter_dataP(1:10:end)),'Color',cmap(3,:),'LineWidth',1)
    xlim([11.3 13.5])
    text(13.4,5,num2str(i),'FontSize',12)
    set(gcf, 'WindowKeyPressFcn', {@selectPeaks}); % 'x' for adding new cycle, 'q' for deleting the last one
    waitfor(gcf);
    
    [~,lks] = findpeaks(idealfilter_data(1:1:end).Data);
    clear x_precision
    for k = 1:numel(pos_manually)
        [~,XI] = min(abs(lks/1e4-pos_manually(k)));
        x_precision(k) = lks(XI);
    end
    
    clear storX storY storZ storQ
    for k = 1:numel(pos_manually)
        storX(k,:) = idealfilter_data.Data(x_precision(k)-249:x_precision(k)+250);
        storY(k,:) = idealfilter_dataP.Data(x_precision(k)-249:x_precision(k)+250);
        storZ(k,:) = A(x_precision(k)-249:x_precision(k)+250);
        storQ(k,:) = B(x_precision(k)-249:x_precision(k)+250);
    end
    if numel(pos_manually) > 0
        selectedCycles{i}.storX = storX;
        selectedCycles{i}.storY = storY;
        selectedCycles{i}.storZ = storZ;
        selectedCycles{i}.storQ = storQ;
        selectedCycles{i}.cellInfo = FileList(i).name(7:10);
        selectedCycles{i}.trial_type = trial_type;
        selectedCycles{i}.x_precision = x_precision;
    else
        selectedCycles{i}.storX = NaN;
        selectedCycles{i}.storY = NaN;
        selectedCycles{i}.storZ = NaN;
        selectedCycles{i}.storQ = NaN;
        selectedCycles{i}.cellInfo = FileList(i).name(7:10);
        selectedCycles{i}.trial_type = trial_type;
    end
end

figure(55)
for k = 1:numel(selectedCycles)
    subplot(ceil(sqrt(numel(selectedCycles))),ceil(sqrt(numel(selectedCycles))),k);
    plot(mean(selectedCycles{k}.storZ)-mean(selectedCycles{k}.storZ(:)),'r');
    hold on; plot(mean(selectedCycles{k}.storQ)-mean(selectedCycles{k}.storQ(:)),'k'); hold off;
    xlabel(selectedCycles{k}.trial_type)
end

figure(33)
for k = 1:numel(selectedCycles)
    clear Sandburg
    try
        for j = 1:size(selectedCycles{k}.storY,1)
            temp = xcorr((selectedCycles{k}.storX(j,:)),(selectedCycles{k}.storY(j,:)));
            temp = temp - min(temp); temp = temp/max(temp);
            Sandburg(j,:) = temp;
        end
        subplot(ceil(sqrt(numel(selectedCycles))),ceil(sqrt(numel(selectedCycles))),k); imagesc(Sandburg)
    end
%     subplot(3,4,k); plot(-499:499,xcorr(mean(selectedCycles{k}.storX(1:13,:)),mean(selectedCycles{k}.storY(1:13,:))))
%     subplot(3,4,k); plot(mean(selectedCycles{k}.storX(1:13,:)),'r'); hold on; plot(mean(selectedCycles{k}.storY(1:13,:)),'k'); hold off
    xlabel(selectedCycles{k}.trial_type)
end

save('selectedCyclesAAAF.mat','selectedCycles','FileList');

% close all
