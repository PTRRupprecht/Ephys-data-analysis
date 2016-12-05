%% go to main folder





%% show electrophysiological recordings
FileList = dir('**.xsg');
try; close 41; close 51; end
offsetV = 0;
offsetI = 0;
suppress_hum = 1;
cmap = lines(numel(FileList));
counter = 1;
VCindex = [];
CCindex = [];
for i = 1:numel(FileList)
    load(FileList(i).name,'-mat');
    A = data.ephys.trace_1;
    

    samplerate = header.ephys.ephys.sampleRate;
    timet = (1:numel(A))/samplerate;%*1000; % ms
    
    
    L_trace = A;
    xdata = [1:size(L_trace)]/samplerate; %in sec
    new_data = timeseries(L_trace,xdata);
    idealfilter_data = idealfilter(new_data,[0 40],'pass');
        
    
    counter = counter + 1;
    figure(41), hold on; plot(idealfilter_data(1:10:end)+offsetV,'Color',cmap(i,:))
    text(timet(end)+1,median(idealfilter_data(1:10:end)+offsetV),num2str(i),'FontSize',12)
    offsetV = offsetV + max(idealfilter_data)-min(idealfilter_data);
    VCindex = [VCindex; i];
end
VCindex = VCindex'
CCindex = CCindex'

