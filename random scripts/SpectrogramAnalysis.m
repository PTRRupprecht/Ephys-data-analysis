


Twenty70 = [21 28 29 30 39 55 63:70 76 86 89 94 ];
Twenty70 = 86;
Twenty00 = [9 17 28 29 38 39 41 42 56 63:70 76 79 82 83 86 89 92 96 98];
Twenty00 = 86;

counter2 = 1;
for jjj = Twenty70
    cd(datasetSingleCells{jjj}.CellID)
    clear S ACF
    %% show electrophysiological recordings
    FileList = dir('**.xsg');
    try; close 42; end
    offsetV = 0;
    cmap = lines(numel(FileList));
    counter = 1;
    setFiles = datasetSingleCells{jjj}.VC70;
    for i = setFiles
        if 1 %strfind(datasetSingleCells{jjj}.VC70odor{counter},'Foo')
            load(FileList(i).name,'-mat');
            A = data.ephys.trace_1;
            if 1 % subtract 50 Hz noise
                window = 10000; % 1 sec
                for kk = 1:numel(A)/window;
                    X = A((1:window) + (kk-1)*window);
                    X = X - mean(X);
                    A((1:window) + (kk-1)*window) = A((1:window) + (kk-1)*window) - repmat(mean( reshape(X,[200 numel(X)/200]),2),[numel(X)/200  1]);
                end
            end
            samplerate = header.ephys.ephys.sampleRate;
            timet = (1:numel(A))/samplerate;%*1000; % ms
            if 0 % subtract 166 Hz noise (?)
                xdata = [1 : size(A,1)] / samplerate; %in sec
                new_data = timeseries(A,xdata);
                idealfilter_data = idealfilter(new_data,[0 100],'pass');
%                 figure(00112); plot(idealfilter_data,'k'); hold on; hold off;
                A = idealfilter_data.Data;
%                 window = 3700; % 1 sec
%                 for kk = 1:numel(A)/window;
%                     X = A((1:window) + (kk-1)*window);
%                     X = X - mean(X);
%                     A((1:window) + (kk-1)*window) = A((1:window) + (kk-1)*window) - repmat(mean( reshape(X,[74 numel(X)/74]),2),[numel(X)/74  1]);
%                 end
            end
%             trace = smooth(A,5);
% 
%             if ~strcmp(header.ephys.ephys.amplifierSettings.Amp_700B_1.mode,'I-Clamp')
%                 figure(42), hold on; plot(timet(1:10:end), trace(1:10:end)+offsetV,'Color',cmap(counter,:))
%                 text(timet(end)+1,median(trace(1:10:end)+offsetV),num2str(i),'FontSize',12)
%                 offsetV = offsetV + max(trace)-min(trace);
%                 xlabel([datasetSingleCells{jjj}.CellID,num2str(jjj)]);
%             end
            window = 15000;
            [S(:,:,counter),F,T] = spectrogram(A(1:30e4),window,window-4000,[3:1:70],1e4);
            ACF(:,counter) = autocorr(A(11e4:14.5e4),1500);
    %         figure(55), subplot(numel(setFiles),1,counter); imagesc(T,F,(mean(abs(S(:,:,counter)),3)))
        end
        counter = counter + 1;
    end
    cd ..
    figure(53),subplot(5,5,counter2); imagesc(T,F,(mean(abs(S(:,:,:)),3))); xlabel(datasetSingleCells{jjj}.CellID,'Interpreter','None');
    figure(545),subplot(5,5,counter2); hold on; plot((1:1501)/10,conv2(ACF(:,:),fspecial('gaussian',[10 1],10),'same'),'r'); xlabel(datasetSingleCells{jjj}.CellID,'Interpreter','None');
    axis([0 150 0.1 1]);
    counter2 = counter2 + 1;
end