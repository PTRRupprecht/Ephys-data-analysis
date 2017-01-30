% s  = short (10:20)
% sm  = shortmedium (10:25)
% ss = very short (7.5:15)
% l  = long (10:40)
% ll = longdelayed (30:60)

% Foo = food, TDC = TDCA, Hi

index = 151 %108
datasetSingleCells{index}.CellID = '161201_AAAD';
datasetSingleCells{index}.pos = [1247 425 79];
datasetSingleCells{index}.posum = [201 -147 60];
datasetSingleCells{index}.VC70 = [1 2 3];
datasetSingleCells{index}.VC70odor = {'Trp' 'Arg' 'Foo'};
datasetSingleCells{index}.PN70III = [];
datasetSingleCells{index}.VC0 = [4:7]; 
datasetSingleCells{index}.VC0odor = {'Foo' 'Foo' 'Foo' 'Foo'};
datasetSingleCells{index}.VCtest = [35 100 50 -5 20];% start sec, ISI msec, duration msec, stim mV, number stim  
datasetSingleCells{index}.CCtest = [];
datasetSingleCells{index}.CCstim = [];
datasetSingleCells{index}.CCodor = {};
datasetSingleCells{index}.possiblyBadTrials = [];
datasetSingleCells{index}.comments = { 'VC time test pulse probing during odor responses, therefore discard many trials (8-10).'};
datasetSingleCells{index}.odors = {'Trp' 'Arg' 'Foo'};
datasetSingleCells{index}.odorLine = [1 2 3];


[ numel(datasetSingleCells{index}.VC70), numel(datasetSingleCells{index}.VC70odor);
 numel(datasetSingleCells{index}.VC0), numel(datasetSingleCells{index}.VC0odor);
numel(datasetSingleCells{index}.CCstim), numel(datasetSingleCells{index}.CCodor)]

datasetSingleCells{index}.pos - datasetSingleCells{index-1}.pos

clear date
cd ..
save(strcat('dataset',date,'_',num2str(index),'.mat'),'datasetSingleCells')

