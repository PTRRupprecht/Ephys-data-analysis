
function [onsetVC70, onsetVC00,odorsVC,dateID] = extractVC(IX,discard,noise50)

clear onsetVC00 onsetVC70

tempPWD = pwd;
cd('C:\Data\rupppete\PhD\electrophysiology2016\SingleCells');
datasetList = dir('dataset*.mat');
load(datasetList(1).name);
dateID = datasetSingleCells{IX}.CellID(1:6);
datasetSingleCells{IX}
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
    for jj = 1:2
        if jj == 1
            choice = trials70;
            indizes = indizes70;
        else
            choice = trials00;
            indizes = indizes00;
        end
        for ii = 1:numel(choice)
            if ~ismember(choice(ii),discard)
                load(traceList(choice(ii)).name,'-mat');
                A = data.ephys.trace_1;
                % remove slow component (exponential)
                if jj == 2
                    ft = fittype( 'a*exp(-x/(b*10000))+c', 'independent', 'x', 'dependent', 'y' );
                    opts = fitoptions( ft ); opts.Display = 'Off';
                    opts.Lower = [-Inf 5 -Inf];
                    opts.StartPoint = [-500 14 -500];
                    opts.Upper = [500 25 500];
                    % Fit model to data.
                    AX = moving(A(1:20:end),200,'median');
                    
                    [fitresult, gof] = fit( (101:numel(A)/20-100)'*20, AX(101:end-100), ft, opts );
                    if abs(fitresult.a) > abs(0.1*fitresult.c)
                        A = (A - fitresult(1:numel(A))) + fitresult.c;
                    end
                end
                if noise50 == 1
                    window = 10000; % 1 sec
                    for kkk = 1:numel(A)/window;
                        X = A((1:window) + (kkk-1)*window);
                        X = X - mean(X);
                        A((1:window) + (kkk-1)*window) = A((1:window) + (kkk-1)*window) - repmat(mean( reshape(X,[200 numel(X)/200]),2),[numel(X)/200  1]);
                    end
                end
                % account for long trials, starting at 30 sec ('llArg', e.g.)
                if jj == 1
                    if strcmp(datasetSingleCells{IX}.VC70odor{indizes(ii)}(1:2),'ll')
                        A = A((2e5+1):end);
                    elseif strcmp(datasetSingleCells{IX}.VC70odor{indizes(ii)}(1:2),'ss')
                        A = [mean(A(1:100))*ones(2.5e4,1);A];
                    end
                else
                    if strcmp(datasetSingleCells{IX}.VC0odor{indizes(ii)}(1:2),'ll')
                        A = A((2e5+1):end);
                    elseif strcmp(datasetSingleCells{IX}.VC0odor{indizes(ii)}(1:2),'ss')
                        A = [mean(A(1:100))*ones(2.5e4,1);A];
                    end
                end
                A(end+1:60e4) = NaN;
                A = A(1:6e5);
                if jj == 1
                    onsetVC70(line_kk,ii,:) = A;
                else
                    onsetVC00(line_kk,ii,:) = A;
                end
            else
                if jj == 1
                    onsetVC70(line_kk,ii,:) = ones(6e5,1)*NaN;
                else
                    onsetVC00(line_kk,ii,:) = ones(6e5,1)*NaN;
                end
            end
        end
    end
end

if ~exist('onsetVC00')
    onsetVC00 = [];
end