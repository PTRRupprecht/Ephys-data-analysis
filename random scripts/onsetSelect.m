
function [onsetVC70, onsetVC00] = onsetSelect(varargin)

if numel(varargin) > 2
    clear keyword
    keyword = varargin{2}.Key;
    gy = varargin{4};
    gx = varargin{3};
    if strcmp(keyword,'x')
        [yyy,xxx]=ginput(1);
        xx = round(xxx); yy = round(yyy);
        distanceVector = sum(([gy;gx]- repmat([xx;yy],[1 size(gy,2)])).^2);
        [~,IX] = min(distanceVector);
        sqrt(min(distanceVector))
    end
else
    IX = varargin{1};
    try; discard = varargin{2}; catch; discard = []; end
end

tempPWD = pwd;
cd('C:\Data\rupppete\PhD\electrophysiology2016\SingleCells');
datasetList = dir('dataset*.mat');
load(datasetList(1).name);

datasetSingleCells{IX}

cd(datasetSingleCells{IX}.CellID) %#ok<USENS>

%% find unique odors for the respective cell (VC)
temp = datasetSingleCells{IX}.VC70odor; for k = 1:numel(temp); temp{k} = temp{k}(2:end); end
for k = 1:numel(temp); if ~ismember(temp{k}(1),'A':'Z');  temp{k} = temp{k}(2:end); end; end
temp2 = datasetSingleCells{IX}.VC0odor; for k = 1:numel(temp2); temp2{k} = temp2{k}(2:end); end
for k = 1:numel(temp2); if ~ismember(temp2{k}(1),'A':'Z');  temp2{k} = temp2{k}(2:end); end; end
odorsVC = unique([temp, temp2]);
%% find unique odors for the respective cell (CC)
temp = datasetSingleCells{IX}.CCodor; for k = 1:numel(temp); temp{k} = temp{k}(2:end); end
for k = 1:numel(temp); if ~ismember(temp{k}(1),'A':'Z');  temp{k} = temp{k}(2:end); end; end
for k = 1:numel(odorsVC); if ~ismember(odorsVC{k}(1),'A':'Z');  odorsVC{k} = odorsVC{k}(2:end); end; end
traceList = dir('*.xsg');

%% plot voltage clamp current traces for different odors
cmap = distinguishable_colors(numel(odorsVC));
offsetV = 0; minVal = 0; maxVal = 0;
minAvg = 0;
maxAvg = 0;


windowmax = 27e4;



A_avg2p = zeros(windowmax,1);
A_avg2m = zeros(windowmax,1);
for kk = 1:numel(odorsVC)
    % check delay
    load('C:\Data\rupppete\PhD\electrophysiology2016\tuningAndTiming\delays04-Aug-2016.mat');
    ff = strfind(datasetSingleCells{IX}.odors,odorsVC{kk});
    delay_index = datasetSingleCells{IX}.odorLine(find(~cellfun(@isempty,ff)));
    odor_delay = 0;% delays(delay_index);
    % find relevant trials
    AA = strfind(datasetSingleCells{IX}.VC70odor,odorsVC{kk});
    indizes = find(~cellfun(@isempty,AA));
    trials70 = datasetSingleCells{IX}.VC70(indizes); %#ok<FNDSB>
    AA = strfind(datasetSingleCells{IX}.VC0odor,odorsVC{kk});
    indizes = find(~cellfun(@isempty,AA));
    trials00 = datasetSingleCells{IX}.VC0(indizes); %#ok<FNDSB>
    figure(701);
    for jj = 1:2
        if jj == 1
            choice = trials70;
        else
            choice = trials00;
        end
        A_avg = zeros(windowmax,1);
        avg_counter = 0;
        for ii = 1:numel(choice)
            if ~ismember(choice(ii),discard)
                load(traceList(choice(ii)).name,'-mat');
                A = data.ephys.trace_1;
                % remove slow component (exponential)
                ft = fittype( 'a*exp(-x/(b*10000))+c', 'independent', 'x', 'dependent', 'y' );
                opts = fitoptions( ft ); opts.Display = 'Off';
                opts.Lower = [-Inf 5 -Inf];
                opts.StartPoint = [-500 14 -500];
                opts.Upper = [500 25 500];
                % Fit model to data.
                [fitresult, gof] = fit( (1:numel(A)/20)'*20, A(1:20:end), ft, opts );
                if abs(fitresult.a) > abs(0.1*fitresult.c)
                    A = (A - fitresult(1:numel(A))) + fitresult.c;
                end


                window = 10000; % 1 sec
                for kkk = 1:numel(A)/window;
                    X = A((1:window) + (kkk-1)*window);
                    X = X - mean(X);
                    A((1:window) + (kkk-1)*window) = A((1:window) + (kkk-1)*window) - repmat(mean( reshape(X,[200 numel(X)/200]),2),[numel(X)/200  1]);
                end
                A = circshift(A,-odor_delay*10);
                samplerate = header.ephys.ephys.sampleRate;
                timet = (1:numel(A))/samplerate;%*1000; % ms
                trace = smooth(A,10);
                figure(701);
                plot(timet(1:10:end), trace(1:10:end)+offsetV,'Color',cmap(kk,:));
                hold on;
                text(timet(end)+1,median(trace(1:10:end)+offsetV),strcat(num2str(choice(ii)),32,odorsVC{kk}),'FontSize',12)
                minVal = min(min(trace(1:10:end)+offsetV),minVal);
                maxVal = max(max(trace(1:10:end)+offsetV),maxVal);

                A_avg = A_avg + A((1:windowmax) + 0e4);
                offsetV = offsetV + max(trace)-min(trace);
                avg_counter = avg_counter + 1;
            end
        end
        drawnow;
        A_avg = A_avg/avg_counter;
        if jj == 1
            A_avg2p = A_avg2p + A_avg;
        else
            A_avg2m = A_avg2m + A_avg;
        end
        minAvg = min(0,min(minAvg,min(A_avg)));
        maxAvg = max(0,max(maxAvg,max(A_avg)));
        
%         figure, plot(cumsum(Y-median(Y)))
        figure(92); plot(timet(1:10:windowmax), smooth(A_avg(1:10:end),10) - mean(A_avg(:)) +(jj-1)*30 - 15,'Color',cmap(kk,:));
        Y = smooth(A_avg(1:1:end),1);
        if jj == 1;
            figure(91); plot(timet(1:1:windowmax), cumsum(Y-mean(Y(1:1e5))),'Color',cmap(kk,:));
            xlim([7.5 16]); 
            [x,y] = ginput(3);
            onsetVC70(kk,:) = x;
            mean(onsetVC70(kk,:),2)
        else
            figure(91); plot(timet(1:1:windowmax), cumsum(Y-mean(Y(1:1e5))),'Color',cmap(kk,:),'LineStyle',':');
            xlim([7.5 16]); 
            [x,y] = ginput(3);
            onsetVC00(kk,:) = x;
            mean(onsetVC00(kk,:),2)
        end    
    end
    offsetV = offsetV + 120;
    figure(701); axis([0 58 minVal-30 maxVal+30]);
    xlabel('time [sec]'); box off; ylabel('current [pA]');
end
figure(92); axis([9 16 minAvg maxAvg]); hold off;
figure(91); xlim([9 16]); hold off;
figure(701); hold off;

datasetSingleCells{IX}.comments{1}

end