%% go to main folder





%% show electrophysiological recordings
FileList = dir('*AAAG*.xsg');
try; close 41; close 51; end
offsetV = 0; offsetV2 = 0;
offsetI = 0;
suppress_hum = 1;
cmap = lines(numel(FileList));
counter = 1;
VCindex = [];
CCindex = [];
counter = 0;
for i = 1:numel(FileList)
    load(FileList(i).name,'-mat');
    samplerate = header.ephys.ephys.sampleRate;
    A = data.ephys.trace_2;
    timet = (1:numel(A))/samplerate;%*1000; % ms

    L_trace = A;
    xdata = [1:size(L_trace)]/samplerate; %in sec
    new_data = timeseries(L_trace,xdata);
    idealfilter_data = idealfilter(new_data,[20 40],'pass');
%     idealfilter_data2 = idealfilter(new_data,[0 60],'pass');
%     idealfilter_data3 = idealfilter(new_data,[0 100],'pass');
%     idealfilter_data4 = idealfilter(new_data,[0 200],'pass');
    
    B = data.ephys.trace_1;
    P_trace = B;
    xdata = [1:size(P_trace)]/samplerate; %in sec
    new_data = timeseries(P_trace,xdata);
%     idealfilter_dataP = idealfilter(new_data,[48 52],'notch');
    idealfilter_dataP = idealfilter(new_data,[0 48],'pass');
    
    
    figure(25), subplot(2,1,1); hold on;% plot(idealfilter_dataP(1:10:end)-mean(idealfilter_dataP(1:10:end)),'Color',cmap(i,:))
    
    plot(idealfilter_data(1:10:end)-mean(idealfilter_data(1:10:end))+counter*10,'Color',cmap(i,:))
%      xlim([65 80])
     figure(25), subplot(2,1,2); hold on;% plot(idealfilter_dataP(1:10:end)-mean(idealfilter_dataP(1:10:end)),'Color',cmap(i,:))
    
    plot(idealfilter_dataP(1:10:end)-mean(idealfilter_dataP(1:10:end))+counter*25,'Color',cmap(i,:))
%     plot(idealfilter_data2(1:10:end)-mean(idealfilter_data2(1:10:end)),'Color',cmap(i+2,:))
%     plot(idealfilter_data3(1:10:end)-mean(idealfilter_data3(1:10:end)),'Color',cmap(i+3,:))
%     plot(idealfilter_data4(1:10:end)-mean(idealfilter_data4(1:10:end)),'Color',cmap(i+4,:))
    
%     akZoom('all')
%     xlim([65 80])
    text(timet(end)+1,0,num2str(i),'FontSize',12)
%     plot(idealfilter_data(1:10:end)-mean(idealfilter_data(1:10:end)),'Color',cmap(i+1,:))
    counter = counter + 1;
% X = smooth(abs(idealfilter_data.Data(600000:900000)),1000);
% Y = smooth(abs(idealfilter_dataP.Data(600000:900000)),1000);
% figure(77), hold on; plot(X(1:1000:end),Y(1:1000:end),'.','Color',cmap(i,:))
% Y = smooth(idealfilter_dataP.Data(600000:900000),100)-  smooth(idealfilter_dataP.Data(600000:900000),20000);
% Z = xcorr(X,Y);
% 
% indx = find(smooth(Y,100)<+5);
% figure, plot(smooth(abs(X(indx)),1000),smooth(Y(indx),1000),'-')
    
    
%     offsetV = offsetV + max(idealfilter_dataP)-min(idealfilter_data); hold off;
    
%     figure(42), hold on; plot(idealfilter_data(1:10:end)+offsetV2,'Color',cmap(i,:))
%     text(timet(end)+1,median(idealfilter_data(1:10:end)+offsetV2),num2str(i),'FontSize',12)
%     offsetV2 = offsetV2 + max(idealfilter_data)-min(idealfilter_data);
%     VCindex = [VCindex; i];
end
VCindex = VCindex'
CCindex = CCindex'


% close all
