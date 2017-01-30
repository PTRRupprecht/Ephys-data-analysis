
function [trialsVC70, trialsVC00] = extractTrialsVC(IX)

clear trialsVC70 trialsVC00

tempPWD = pwd;
cd('C:\Data\rupppete\PhD\electrophysiology2016\SingleCells');
datasetList = dir('dataset*.mat');
load(datasetList(1).name);
cd(datasetSingleCells{IX}.CellID)


%% find unique odors for the respective cell (VC)
temp = datasetSingleCells{IX}.VC70odor; for k = 1:numel(temp); temp{k} = temp{k}(2:end); end
for k = 1:numel(temp); if ~ismember(temp{k}(1),'A':'Z');  temp{k} = temp{k}(2:end); end; end
temp2 = datasetSingleCells{IX}.VC0odor; for k = 1:numel(temp2); temp2{k} = temp2{k}(2:end); end
for k = 1:numel(temp2); if ~ismember(temp2{k}(1),'A':'Z');  temp2{k} = temp2{k}(2:end); end; end
odorsVC = unique([temp, temp2]);
traceList = dir('*.xsg');

%% plot voltage clamp current traces for different odors

for kk = 1:numel(odorsVC)
    % check delay
    ff = strfind(datasetSingleCells{IX}.odors,odorsVC{kk});
    line_kk = datasetSingleCells{IX}.odorLine(find(~cellfun(@isempty,ff)));
    
    % find relevant trials
    AA = strfind(datasetSingleCells{IX}.VC70odor,odorsVC{line_kk});
    indizes70 = find(~cellfun(@isempty,AA));
    trials70 = datasetSingleCells{IX}.VC70(indizes70); %#ok<FNDSB>
    AA = strfind(datasetSingleCells{IX}.VC0odor,odorsVC{line_kk});
    indizes00 = find(~cellfun(@isempty,AA));
    trials00 = datasetSingleCells{IX}.VC0(indizes00); %#ok<FNDSB>
    
    trialsVC70(line_kk,1:numel(trials70)) = trials70;
    trialsVC00(line_kk,1:numel(trials00)) = trials00;
end

trialsVC00(trialsVC00 == 0) = NaN;
trialsVC70(trialsVC70 == 0) = NaN;

cd(tempPWD)
