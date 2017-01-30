
function [odorsVCII,odorLine] = extractOdorsAndLines(IX)

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
    clear odorsVCII
    for kk = 1:numel(odorsVC)
        ff = strfind(datasetSingleCells{IX}.odors,odorsVC{kk});
        line_kk = datasetSingleCells{IX}.odorLine(find(~cellfun(@isempty,ff)));
        odorsVCII{kk} = odorsVC{line_kk};
    end
    for m = 1:numel(odorsVC)
        odorLine(m) = find( strcmp(odorsVC{m},datasetSingleCells{IX}.odors));
    end
    cd(tempPWD)
end