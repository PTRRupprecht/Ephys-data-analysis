%% go to main folder





%% show electrophysiological recordings
FileList = dir('*AAAH*.xsg');
try; close 41; close 51; end
offsetV = 0; offsetV2 = 0;
offsetI = 0;
suppress_hum = 1;
cmap = lines(numel(FileList));
counter = 1;
VCindex = [];
CCindex = [];
for i = 1:numel(FileList)
    load(FileList(i).name,'-mat');
    samplerate = header.ephys.ephys.sampleRate;
    timet = (1:numel(A))/samplerate;%*1000; % ms

    A = data.ephys.trace_2;
    L_trace = A;
    xdata = [1:size(L_trace)]/samplerate; %in sec
    new_data = timeseries(L_trace,xdata);
    idealfilter_data = idealfilter(new_data,[0 40],'pass');
%     idealfilter_data2 = idealfilter(new_data,[0 60],'pass');
%     idealfilter_data3 = idealfilter(new_data,[0 100],'pass');
%     idealfilter_data4 = idealfilter(new_data,[0 200],'pass');
    
    B = data.ephys.trace_1;
    P_trace = B;
    xdata = [1:size(P_trace)]/samplerate; %in sec
    new_data = timeseries(P_trace,xdata);
%     idealfilter_dataP = idealfilter(new_data,[48 52],'notch');
    idealfilter_dataP = idealfilter(new_data,[0 48],'pass');
    
    
    counter = counter + 1;
    figure(24), hold on;% plot(idealfilter_dataP(1:10:end)-mean(idealfilter_dataP(1:10:end)),'Color',cmap(i,:))
    
    plot(idealfilter_data(1:10:end)-mean(idealfilter_data(1:10:end))+i*10,'Color',cmap(i,:))
    figure(25), hold on;% plot(idealfilter_dataP(1:10:end)-mean(idealfilter_dataP(1:10:end)),'Color',cmap(i,:))
    
    plot(idealfilter_dataP(1:10:end)-mean(idealfilter_dataP(1:10:end))+i*50,'Color',cmap(i,:))
%     plot(idealfilter_data2(1:10:end)-mean(idealfilter_data2(1:10:end)),'Color',cmap(i+2,:))
%     plot(idealfilter_data3(1:10:end)-mean(idealfilter_data3(1:10:end)),'Color',cmap(i+3,:))
%     plot(idealfilter_data4(1:10:end)-mean(idealfilter_data4(1:10:end)),'Color',cmap(i+4,:))
    
    
    text(timet(end)+1,0,num2str(i),'FontSize',12)
%     plot(idealfilter_data(1:10:end)-mean(idealfilter_data(1:10:end)),'Color',cmap(i+1,:))
    
    
    
%     offsetV = offsetV + max(idealfilter_dataP)-min(idealfilter_data); hold off;
    
%     figure(42), hold on; plot(idealfilter_data(1:10:end)+offsetV2,'Color',cmap(i,:))
%     text(timet(end)+1,median(idealfilter_data(1:10:end)+offsetV2),num2str(i),'FontSize',12)
%     offsetV2 = offsetV2 + max(idealfilter_data)-min(idealfilter_data);
%     VCindex = [VCindex; i];
end
VCindex = VCindex'
CCindex = CCindex'


% close all
