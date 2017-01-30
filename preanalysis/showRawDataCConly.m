
function [VCtraces] = showRawDataCConly(varargin)

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
    discard = [];
else
    IX = varargin{1};
    try; discard = varargin{2}; catch; discard = []; end
end

tempPWD = pwd;
cd('C:\Data\rupppete\PhD\electrophysiology2016\SingleCells');
datasetList = dir('dataset*.mat');
load(datasetList(1).name);

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

numTimes = 27e4;

cmap = distinguishable_colors(3);
%% plot current clamp voltage traces for different odors
offsetV = 0; minVal = 0; maxVal = 0;
if ~isempty(datasetSingleCells{IX}.CCodor)
    for kk = 1:numel(odorsVC)
        % check delay
        load('C:\Data\rupppete\PhD\electrophysiology2016\tuningAndTiming\delays04-Aug-2016.mat');
        ff = strfind(datasetSingleCells{IX}.odors,odorsVC{kk});
        delay_index = datasetSingleCells{IX}.odorLine(find(~cellfun(@isempty,ff)));
        odor_delay = delays(delay_index);
        % find relevant trials
        AA = strfind(datasetSingleCells{IX}.CCodor,odorsVC{kk});
        indizes = find(~cellfun(@isempty,AA));
        trialsCC = datasetSingleCells{IX}.CCstim(indizes);  %#ok<FNDSB>
        figure(702);
        choice = trialsCC;
        for ii = 1:numel(trialsCC)
            load(traceList(choice(ii)).name,'-mat');
            A = data.ephys.trace_1;
            window = 10000; % 1 sec
            for kkk = 1:numel(A)/window;
                X = A((1:window) + (kkk-1)*window);
                X = X - mean(X);
                A((1:window) + (kkk-1)*window) = A((1:window) + (kkk-1)*window) - repmat(mean( reshape(X,[200 numel(X)/200]),2),[numel(X)/200  1]);
            end
%             A = circshift(A,-odor_delay*10);
            samplerate = header.ephys.ephys.sampleRate;
            timet = (1:numel(A))/samplerate;%*1000; % ms
            trace = smooth(A,1);
            plot(timet(1:1:end), trace(1:1:end)+offsetV,'Color',cmap(kk,:));
            hold on;
            text(timet(end)+1,median(trace(1:1:end)+offsetV),strcat(num2str(choice(ii)),32,odorsVC{kk}),'FontSize',12)
            minVal = min(min(trace(1:1:end)+offsetV),minVal);
            maxVal = max(max(trace(1:1:end)+offsetV),maxVal);

            offsetV = offsetV + max(trace)-min(trace);
        end
        axis([0 58 minVal-10 maxVal+10]);
        xlabel('time [sec]'); box off; ylabel('voltage [mV]');
    end
    hold off;
else
    figure(702); plot(0,0);
end

datasetSingleCells{IX}

end