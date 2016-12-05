
function showRawDataOnlyVC(varargin)

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
traceList = dir('*.xsg');


%% plot voltage clamp current traces for different odors
cmap = distinguishable_colors(numel(odorsVC));
offsetV = 0; minVal = 0; maxVal = 0;
for kk = 1:numel(odorsVC)
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
        for ii = 1:numel(choice)
            if ~ismember(choice(ii),discard)
                load(traceList(choice(ii)).name,'-mat');
                A = data.ephys.trace_1;
                samplerate = header.ephys.ephys.sampleRate;
                timet = (1:numel(A))/samplerate;%*1000; % ms
                trace = smooth(A,10);
                figure(701);
                plot(timet(1:10:end), trace(1:10:end)+offsetV,'Color',cmap(kk,:));
                hold on;
                text(timet(end)+1,median(trace(1:10:end)+offsetV),strcat(num2str(choice(ii)),32,odorsVC{kk}),'FontSize',12)
                minVal = min(min(trace(1:10:end)+offsetV),minVal);
                maxVal = max(max(trace(1:10:end)+offsetV),maxVal);

                offsetV = offsetV + max(trace)-min(trace);
            end
        end
    end
    offsetV = offsetV + 120;
    figure(701); axis([0 58 minVal-30 maxVal+30]);
    xlabel('time [sec]'); box off; ylabel('current [pA]');
end
figure(701); hold off;

datasetSingleCells{IX}
try
datasetSingleCells{IX}.comments{1}
end
end