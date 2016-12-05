

%% Manual selection of bad trials
for IX = 148:150% 1:numel(datasetSingleCells)
    
    IX
    clear onsetVC00 onsetVC70
    showRawDataOnlyVC(IX);

    answer = inputdlg({'Not good for transients:','Not good elsewhere:','Comments'},'UI / Trial selection',1);

    discard = str2num(answer{1});
    discard_later = str2num(answer{2});
    
    [onsetVC70, onsetVC00,odorsVC,dateID] = extractVC(IX,discard,1);
    
    X_all{IX}.onsetVC70 = onsetVC70;
    X_all{IX}.onsetVC00 = onsetVC00;
    X_all{IX}.dateID = dateID;
    X_all{IX}.discard_later = discard_later;
    X_all{IX}.comments = str2num(answer{3});
    % discarded trials: look for NaNs ... (repair this issue later!)

    [odorsVC,odorLine] = extractOdorsAndLines(IX);
    sTraces{IX}.odorsVCII = odorsVC;
    sTraces{IX}.odorLine = odorLine;
    
end
sTraces = X_all;
save('datasetSelectedTracesCorrected.mat','sTraces','-v7.3')

%% Manual selection of response onset for each odor/fish combination
% full automatisation does not make a lot of sense. Sometimes the onset is
% simply not clear and onset should be rather used from fish taken at
% adjacent days. Also, flow changes due to pump wheel speed changes can
% induce responses. Those, however, occur at the same time for all
% odorants, and are to be ignored in this odor onset timing analysis,
% although the exact determination of odor onset may be difficult if it
% overlaps with the response to the flow change response.

% neurons (dateXnew) and fish (dateUnique)
counter = 1; nennersum = 0; totalsum = 0;
cmap = distinguishable_colors(3);
for j = 1:numel(datasetSingleCells)
    dateXnew(j) = str2double(datasetSingleCells{j}.CellID(1:6));
end
dateUnique = unique(dateXnew); dateUnique = sort(dateUnique);

% average across all trials and neurons for each odor/fish combination for later use in GUI
for k = 33 % 1:numel(dateUnique)
    cells2do = find(dateXnew == dateUnique(k));
    for jj = 1:3; super70{k,jj} = zeros(60e4,1); super00{k,jj} = zeros(60e4,1); end
    for cc = 1:numel(cells2do)
        for jj = 1:size(sTraces{cells2do(cc)}.onsetVC70,1)
            try
                if ~isnan(nanmean(squeeze(nanmean(sTraces{cells2do(cc)}.onsetVC70(jj,:,:),2))))
                    super70{k,jj} = super70{k,jj} + squeeze(nanmean(sTraces{cells2do(cc)}.onsetVC70(jj,:,:),2));
                end
            end
            try
                if ~isnan(nanmean(squeeze(nanmean(sTraces{cells2do(cc)}.onsetVC00(jj,:,:),2))))
                    super00{k,jj} = super00{k,jj} + squeeze(nanmean(sTraces{cells2do(cc)}.onsetVC00(jj,:,:),2));
                end
            end
        end
    end
end

% graphical representation manual choice of onset (use excel file
% (Response ... .xlsx) to write down values)
k = 33 % fish number to analyze
cells2do = find(dateXnew == dateUnique(k));
cells2do % indicates number of cells for this fish
cmap = distinguishable_colors(3);
for jj = 1:3
    figure(4);
    plot((7e4:18e4)/1e4,smooth(super70{k,jj}(7e4:18e4,1),60) - max(smooth(super70{k,jj}(1:26e4,1),40)),'Color',cmap(jj,:)); hold on;
    plot((7e4:18e4)/1e4,smooth(super00{k,jj}(7e4:18e4,1),60) - min(smooth(super00{k,jj}(1:26e4,1),40)),'Color',cmap(jj,:)); hold on;
end
xlabel(num2str(dateUnique(k))); hold off;
figure(1); j = 1;
plot((7e4:18e4)/1e4,smooth(super70{k,j}(7e4:18e4,1),60) - max(smooth(super70{k,j}(1:26e4,1),40)),'Color',cmap(j,:)); hold on;
plot((7e4:18e4)/1e4,smooth(super00{k,j}(7e4:18e4,1),60) - min(smooth(super00{k,j}(1:26e4,1),40)),'Color',cmap(j,:)); hold off;
figure(2); j = 2;
plot((7e4:18e4)/1e4,smooth(super70{k,j}(7e4:18e4,1),60) - max(smooth(super70{k,j}(1:26e4,1),40)),'Color',cmap(j,:)); hold on;
plot((7e4:18e4)/1e4,smooth(super00{k,j}(7e4:18e4,1),60) - min(smooth(super00{k,j}(1:26e4,1),40)),'Color',cmap(j,:)); hold off;
figure(3); j = 3;
plot((7e4:18e4)/1e4,smooth(super70{k,j}(7e4:18e4,1),60) - max(smooth(super70{k,j}(1:26e4,1),40)),'Color',cmap(j,:)); hold on;
plot((7e4:18e4)/1e4,smooth(super00{k,j}(7e4:18e4,1),60) - min(smooth(super00{k,j}(1:26e4,1),40)),'Color',cmap(j,:)); hold off;

% check out single cells (example)
showRawDataOnlyVC(87)

% first import ('ResponseOnsetPerFish' (Excel) as matrix)
% then plot delays per fish
for k = 1:numel(dateUnique)
    cells2do = find(dateXnew == dateUnique(k));
    [~,sIX] = sort(sTraces{cells2do(1)}.odorLine); % overkill - sIX is always the same as the odorLine variable ...
    if k < 11 % manual correction for different line delays for those experiments !!!
        sIX = [3 1 2];
    end
    ResponseTime(k,1) = ResponseOnsetPerFish(k,sIX(1)+2);
    if numel(sIX)>1; ResponseTime(k,2) = ResponseOnsetPerFish(k,sIX(2)+2); end
    if numel(sIX)>2; ResponseTime(k,3) = ResponseOnsetPerFish(k,sIX(3)+2); end
end
ResponseTime(ResponseTime==0) = NaN;

% fill in and interpolate the missing response time values
for p = 1:size(ResponseTime,1)
    onsetX2(p,:) = ResponseTime(p,:) - nanmean(ResponseTime(p,:));
end
onset1 = nanmean(onsetX2(:,1));
onset2 = nanmean(onsetX2(:,2));
onset3 = nanmean(onsetX2(:,3));

for p = 1:size(ResponseTime,1)
    onsetZ(p,:) = nanmean(ResponseTime(p,:) - [onset1 onset2 onset3]) + [onset1 onset2 onset3];
end
onsetZ(5,:) = NaN;
for j = 1:3
    good_indizes = find(~isnan(onsetZ(:,j)));
    onsetZ(:,j) = interp1(good_indizes,onsetZ(good_indizes,j),1:size(onsetZ,1),'Nearest','extrap');
end
ResponseTimeAll = onsetZ;

% save response times in the common structure
for IX = 1:numel(datasetSingleCells)
    k = find(dateXnew(IX) == dateUnique);
    cells2do = find(dateXnew == dateUnique(k));
    [~,sIX] = sort(sTraces{cells2do(1)}.odorLine); % overkill, again
    if k < 11 % manual correction for different line delays for those experiments !!!
        sIX = [2 3 1];
    end
    for j = 1:3
        XIs = find(sIX == j);
        ssTraces{IX}.delay(XIs) = ResponseTimeAll(k,j);
    end
end

%% retrospectively find out trial numbers of all timetraces
for k =148:150 %1:numel(datasetSingleCells)
    [trialsVC70, trialsVC00] = extractTrialsVC(k);
    ssTraces{k}.trialsVC70 = trialsVC70;
    ssTraces{k}.trialsVC00 = trialsVC00;
end


%% select neurons to be discarded for tuning analysis (decreasing access resistance)

for IX = 148:150
    stableRRR{IX} = manuallyDeselectRacc(IX,datasetSingleCells{IX},sTraces{IX},ssTraces{IX});
end

for k = 148:150%1:numel(stableRRR)
    ssTraces{k}.stableResponses = stableRRR{k};
end


save('dataSelected_metaX.mat','ssTraces')

%% select brain region for each neuron

load('C:\Data\rupppete\PhD\electrophysiology2016\AnnotationViewer\StacksUint16.mat');
figure(99); imshow3DofDp(Horizontal,Sagittal,'NeuronList_Annotation',[8387 8557]);

IX = 148;

IX = IX + 1

NeuronList{1}.x = datasetSingleCells{IX}.pos(1);
NeuronList{1}.y = datasetSingleCells{IX}.pos(2);
NeuronList{1}.z = datasetSingleCells{IX}.pos(3);
NeuronList{1}.x_abs = datasetSingleCells{IX}.posum(1);
NeuronList{1}.y_abs = datasetSingleCells{IX}.posum(2);
NeuronList{1}.z_abs = datasetSingleCells{IX}.posum(3);
NeuronList{1}.date = '';
NeuronList{1}.comment = 'nothing';
save('NeuronList_Annotation.mat','NeuronList');


%% import anatomical locations overview as excel sheet >> matrix
% and show them on the sagittal, horizontal and surfacial projections
nb_Dp = nansum(Anatomicallocationsoverview(:,2)./nansum(Anatomicallocationsoverview(:,2:end),2));
nb_DpNTborder = nansum(Anatomicallocationsoverview(:,3)./nansum(Anatomicallocationsoverview(:,2:end),2));
nb_NT = nansum(Anatomicallocationsoverview(:,4)./nansum(Anatomicallocationsoverview(:,2:end),2));
nb_furrow = nansum(Anatomicallocationsoverview(:,5)./nansum(Anatomicallocationsoverview(:,2:end),2));
nb_aDpfurrowborder = nansum(Anatomicallocationsoverview(:,6)./nansum(Anatomicallocationsoverview(:,2:end),2));
nb_aDp = nansum(Anatomicallocationsoverview(:,7)./nansum(Anatomicallocationsoverview(:,2:end),2));
nb_notknown = nansum(Anatomicallocationsoverview(:,8)./nansum(Anatomicallocationsoverview(:,2:end),2));
% sanity check
nb_total = nb_Dp + nb_DpNTborder + nb_NT + nb_furrow + nb_aDpfurrowborder + nb_aDp + nb_notknown


[anatomicalColorsY,anatomicalColorsX] = find(~isnan(Anatomicallocationsoverview(:,2:end)));
NeuronViewerMap(anatomicalColorsX,anatomicalColorsY,51324)





%% real data analysis

numel00 = 0; numel70 = 0;
for cellIX = 1:numel(sTraces)
    for odorIX = 1:size(sTraces{cellIX}.onsetVC70,1)
        numel00 = numel00+size(sTraces{cellIX}.onsetVC00,2);
        numel70 = numel70+size(sTraces{cellIX}.onsetVC70,2);
    end
end

criterion = Anatomicallocationsoverview(:,2) == 1;% | Anatomicallocationsoverview(:,3) == 1; % only Dp
% criterion = any(Anatomicallocationsoverview(:,2:8)' == 1); % any cell

counter70 = 1; counter00 = 1; counterOdor = 1;
trace00all = zeros(numel00,60000); trace70all = zeros(numel70,60000); % all trials, odors, neurons
trace70odor = zeros(400,60000); trace00odor = zeros(400,60000); % average for neuron-odor pairs across trials; all odors, neurons
smoothing = fspecial('gaussian',[10 1],5);
for cellIX = 1:numel(sTraces)
    cellIX
    if criterion(cellIX)
        for odorIX = 1:size(sTraces{cellIX}.onsetVC70,1)
            odorCounter00 = 0;
            for trialIX = 1:size(sTraces{cellIX}.onsetVC00,2)
                trace00 = circshift(squeeze(sTraces{cellIX}.onsetVC00(odorIX,trialIX,:)),[-round(ssTraces{cellIX}.delay(odorIX)*1e4)+10e4 0]);
                if ~isnan(trace00(1)) && trace00(1) && any(str2num(ssTraces{cellIX}.stableResponses{2})==ssTraces{cellIX}.trialsVC00(odorIX,trialIX))
                    trace00 = conv(trace00,smoothing,'same');
                    trace00 = trace00(1:10:end);
                    trace00 = (trace00 - nanmean(trace00(1:0.8e4)))/nanstd(trace00(1:0.8e4));
                    trace00all(counter00,:) = trace00;
                    counter00 = counter00 + 1;
                    odorCounter00 = odorCounter00 + 1;
                end
            end
            trace00odor(counterOdor,:) = nanmean(trace00all(counter00-odorCounter00+1:counter00,:),1);
            odorCounter70 = 0;
            for trialIX = 1:size(sTraces{cellIX}.onsetVC70,2)
                trace70 = circshift(squeeze(sTraces{cellIX}.onsetVC70(odorIX,trialIX,:)),[-round(ssTraces{cellIX}.delay(odorIX)*1e4)+10e4 0]);
                if ~isnan(trace70(1)) && trace70(1) && any(str2num(ssTraces{cellIX}.stableResponses{1})==ssTraces{cellIX}.trialsVC70(odorIX,trialIX))
                    trace70 = conv(trace70,smoothing,'same');
                    trace70 = trace70(1:10:end);
                    trace70 = (trace70 - nanmean(trace70(1:0.8e4)))/nanstd(trace70(1:0.8e4));
                    trace70all(counter70,:) = trace70;
                    counter70 = counter70 + 1;
                    odorCounter70 = odorCounter70 + 1;
                end
            end
            trace70odor(counterOdor,:) = nanmean(trace70all(counter70-odorCounter70+1:counter70,:),1);
            counterOdor = counterOdor + 1;
        end
    end
end


chosenNeuronOdors = find(~isnan(trace00odor(:,1e4)) & trace00odor(:,1e4)~=0 & ~isnan(trace70odor(:,1e4)) & trace70odor(:,1e4)~=0 );

chosen00odor = trace00odor(chosenNeuronOdors,:);
chosen70odor = trace70odor(chosenNeuronOdors,:);

rankingResponse = nanmean(chosen00odor(:,1.0e4:1.15e4),2);
rankingResponse = nanmean(chosen70odor(:,1.3e4:2e4),2);
rankingResponse = nanmean(chosen70odor(:,2.1e4:2.5e4),2);
[~,XI] = sort(rankingResponse,'descend');
figure(1),subplot(1,2,1); imagesc([0 60],[],chosen70odor(XI,:),[-2 2]); xlim([5 20]);
xlabel('time (sec)');
 subplot(1,2,2);imagesc([0 60],[],chosen00odor(XI,:),[-2 2]); xlim([5 20]);
xlabel('time (sec)');

timeT = (1:6e4)/1e3;
figure(3), plot(timeT,mean(chosen00odor),'r'); hold on; plot(timeT,-mean(chosen70odor),'k')
xlabel('time (sec)'); ylabel('Currents (normalized, then averaged)'); xlim([9 20]); hold off

%% quantify area-under-curve after supposed response onset in a 1.5 sec-time window

criterion = Anatomicallocationsoverview(:,2) == 1;% | Anatomicallocationsoverview(:,3) == 1; % only Dp
% criterion = any(Anatomicallocationsoverview(:,2:8)' == 1); % any cell

counter70 = 1; counter00 = 1; counterOdor = 1;
trace00all = zeros(numel00,60000); trace70all = zeros(numel70,60000); % all trials, odors, neurons
trace70odor = zeros(400,60000); trace00odor = zeros(400,60000); % average for neuron-odor pairs across trials; all odors, neurons
smoothing = fspecial('gaussian',[10 1],5);
for cellIX = 1:numel(sTraces)
    cellIX
    if criterion(cellIX)
        for odorIX = 1:size(sTraces{cellIX}.onsetVC70,1)
            odorCounter00 = 0;
            for trialIX = 1:size(sTraces{cellIX}.onsetVC00,2)
                trace00 = circshift(squeeze(sTraces{cellIX}.onsetVC00(odorIX,trialIX,:)),[-round(ssTraces{cellIX}.delay(odorIX)*1e4)+10e4 0]);
                if ~isnan(trace00(1)) && trace00(1) && any(str2num(ssTraces{cellIX}.stableResponses{2})==ssTraces{cellIX}.trialsVC00(odorIX,trialIX))
                    trace00 = conv(trace00,smoothing,'same');
                    trace00 = trace00(1:10:end);
                    trace00 = (trace00 - nanmean(trace00(1:0.8e4)));% /nanstd(trace00(1:0.8e4));
                    trace00all(counter00,:) = trace00;
                    counter00 = counter00 + 1;
                    odorCounter00 = odorCounter00 + 1;
                end
            end
            trace00odor(counterOdor,:) = nanmean(trace00all(counter00-odorCounter00+1:counter00,:),1);
            odorCounter70 = 0;
            for trialIX = 1:size(sTraces{cellIX}.onsetVC70,2)
                trace70 = circshift(squeeze(sTraces{cellIX}.onsetVC70(odorIX,trialIX,:)),[-round(ssTraces{cellIX}.delay(odorIX)*1e4)+10e4 0]);
                if ~isnan(trace70(1)) && trace70(1) && any(str2num(ssTraces{cellIX}.stableResponses{1})==ssTraces{cellIX}.trialsVC70(odorIX,trialIX))
                    trace70 = conv(trace70,smoothing,'same');
                    trace70 = trace70(1:10:end);
                    trace70 = (trace70 - nanmean(trace70(1:0.8e4)));%/nanstd(trace70(1:0.8e4));
                    trace70all(counter70,:) = trace70;
                    counter70 = counter70 + 1;
                    odorCounter70 = odorCounter70 + 1;
                end
            end
            trace70odor(counterOdor,:) = nanmean(trace70all(counter70-odorCounter70+1:counter70,:),1);
            counterOdor = counterOdor + 1;
        end
    end
end

chosenNeuronOdors = find(~isnan(trace00odor(:,1e4)) & trace00odor(:,1e4)~=0 & ~isnan(trace70odor(:,1e4)) & trace70odor(:,1e4)~=0 );
chosen00odor = trace00odor(chosenNeuronOdors,:);
chosen70odor = trace70odor(chosenNeuronOdors,:);

integWindow = 1e4:1.15e4;
AUC_00 = mean(chosen00odor(:,integWindow)')*numel(integWindow)/1e3;
AUC_70 = mean(chosen70odor(:,integWindow)')*numel(integWindow)/1e3;

figure(4); plot(AUC_00,-AUC_70,'.g','Markersize',16)
axis([-5 45 -5 37.5])


%% phase plot for exc. and inhibition (trajectories in EPSC/IPSC space)

cmap = jet(150);
counter = 0;
ch70All = [];%zeros(1501,1);
ch00All = [];%zeros(1501,1);
figure(54), 
for k = 1:size(chosen00odor,1)
    ch00 = smooth(chosen00odor(k,1e4:2e4),150);
    ch70 = -smooth(chosen70odor(k,1e4:2e4),150);
    if max(ch00) > 5 && max(ch70) > 5
        counter = counter + 1;
        ch00 = (ch00 - min(ch00))/(max(ch00) - min(ch00));
        ch00All = [ch00All, ch00];
        ch70 = (ch70 - min(ch70))/(max(ch70) - min(ch70));
        ch70All = [ch70All, ch70];
        for j = 1:50:9900%890:-10:1%1:10:1500
            plot(ch70((j):(j+51)),ch00((j):(j+51)),'Color',cmap(min(150,ceil(j/40)),:)); hold on;
        end
    end
end;
ch70All = mean(ch70All');
ch00All = mean(ch00All');
for j = 1:50:9900%890:-10:1%1:10:1500
    plot(ch70All((j):(j+51)),ch00All((j):(j+51)),'Color',cmap(min(150,ceil(j/40)),:),'LineWidth',5); hold on;
end
hold off



%% tuning a la Thomas (plus: scrambled for control)
% select region ('criterion'), inh/exc ('traceXXall' etc.) and the
% corresponding plotting modality ('order' etc., and subplot panel) and
% possiblythe integration window of interest ('integWindow')


criterion = Anatomicallocationsoverview(:,2) == 1;% | Anatomicallocationsoverview(:,3) == 1; % only Dp
% criterion = any(Anatomicallocationsoverview(:,2:8)' == 1); % any cell

integWindow = 1e4:1.15e4;
% integWindow = 1.3e4:2e4;
counter70 = 1; counter00 = 1; counterOdor = 1;
trace00all = zeros(numel00,1); trace70all = zeros(numel70,1); % all trials, odors, neurons
smoothing = fspecial('gaussian',[10 1],5);
clear trace70all trace70odor trace70cell trace00all trace00odor trace00cell
for cellIX = 1:numel(sTraces)
    cellIX
    if criterion(cellIX)
        for odorIX = 1:size(sTraces{cellIX}.onsetVC70,1)
            odorCounter00 = 0;
            for trialIX = 1:size(sTraces{cellIX}.onsetVC00,2)
                trace00 = circshift(squeeze(sTraces{cellIX}.onsetVC00(odorIX,trialIX,:)),[-round(ssTraces{cellIX}.delay(odorIX)*1e4)+10e4 0]);
                if ~isnan(trace00(1)) && trace00(1) && any(str2num(ssTraces{cellIX}.stableResponses{2})==ssTraces{cellIX}.trialsVC00(odorIX,trialIX))
                    trace00 = conv(trace00,smoothing,'same');
                    trace00 = trace00(1:10:end);
                    trace00 = (trace00 - nanmean(trace00(1:0.8e4)));% /nanstd(trace00(1:0.8e4));
                    
                    
                    trace00all(counter00) = mean(trace00(integWindow));
                    trace00odor(counter00) = odorIX;
                    trace00cell(counter00) = cellIX;
                    counter00 = counter00 + 1;
                    
                end
            end
            for trialIX = 1:size(sTraces{cellIX}.onsetVC70,2)
                trace70 = circshift(squeeze(sTraces{cellIX}.onsetVC70(odorIX,trialIX,:)),[-round(ssTraces{cellIX}.delay(odorIX)*1e4)+10e4 0]);
                if ~isnan(trace70(1)) && trace70(1) && any(str2num(ssTraces{cellIX}.stableResponses{1})==ssTraces{cellIX}.trialsVC70(odorIX,trialIX))
                    trace70 = conv(trace70,smoothing,'same');
                    trace70 = trace70(1:10:end);
                    trace70 = (trace70 - nanmean(trace70(1:0.8e4)));%/nanstd(trace70(1:0.8e4));
                    trace70all(counter70) = mean(trace70(integWindow));
                    trace70odor(counter70) = odorIX;
                    trace70cell(counter70) = cellIX;
                    counter70 = counter70 + 1;
                end
            end
        end
    end
end

traceXXall = trace70all;
traceXXodor = trace70odor;
traceXXcell = trace70cell;
traceXXall = trace00all;
traceXXodor = trace00odor;
traceXXcell = trace00cell;

counterCell = 1;
clear AUCs2 AUCs2multi AUCs2permute AUCspermute
for k = 1:numel(sTraces)
    trials = find(traceXXcell == k);
    if ~isempty(trials)
        AUCs = traceXXall(trials);
        for pp = 1:100
            AUCspermute{pp} = AUCs(randperm(numel(AUCs(:))));
        end
        odors = traceXXodor(trials);
        if numel(unique(odors)) == 3
            for jj = 1:3
                trials2 = find(odors == jj);
                AUCs2(counterCell,jj) = mean(AUCs(trials2));
                for pp = 1:100
                    AUCs2permute{pp}(counterCell,jj) = mean(AUCspermute{pp}(trials2));
                end
                AUCs2multi(counterCell,jj) = numel(trials2);
            end
            counterCell = counterCell + 1;
        end
    end
end
AUCs2 = AUCs2'; AUCs2multi = AUCs2multi';
clear AUCs2multiSort
[AUCs2sort,XI] = sort(AUCs2);
for j = 1:size(XI,2)
    AUCs2multiSort(:,j) = AUCs2multi(XI(:,j),j);
end
clear AUCs2sort2
for pp = 1:100
    AUCs2permute{pp} = AUCs2permute{pp}';
    clear AUCs2multiSort2
    [AUCs2sort2{pp},XI] = sort(AUCs2permute{pp});
    for j = 1:size(XI,2)
        AUCs2multiSort2{pp}(:,j) = AUCs2multi(XI(:,j),j);
    end
end

unscrambled = mean(AUCs2sort');
scrambled = mean(cell2mat(AUCs2sort2)');

order = [3 2 1];
factorSign = +1;
thickness = +4;
order = [1 2 3];
factorSign = -1;
thickness = +2;

figure(7), subplot(1,4,1);
plot(factorSign*unscrambled(order),'LineWidth',thickness); hold on; plot(factorSign*scrambled(order),'r','LineWidth',thickness);
xlim([0.8 3.2]); %ylim([2 16])



%% variability explained by a) tuning b) adaptation c) stimulus strength


criterion = Anatomicallocationsoverview(:,2) == 1;% | Anatomicallocationsoverview(:,3) == 1; % only Dp
% criterion = any(Anatomicallocationsoverview(:,2:8)' == 1); % any cell

integWindow = 1e4:1.15e4;
% integWindow = 1.4e4:2e4;
counter70 = 1; counter00 = 1; counterOdor = 1;
trace00all = zeros(numel00,1); trace70all = zeros(numel70,1); % all trials, odors, neurons
smoothing = fspecial('gaussian',[10 1],5);
clear trace70all trace70odor trace70cell trace00all trace00odor trace00cell trace70trial trace00trial trace00odorid trace70odorid
for cellIX = 1:numel(sTraces)
    cellIX
    if criterion(cellIX)
        for odorIX = 1:size(sTraces{cellIX}.onsetVC70,1)
            odorCounter00 = 0;
            for trialIX = 1:size(sTraces{cellIX}.onsetVC00,2)
                trace00 = circshift(squeeze(sTraces{cellIX}.onsetVC00(odorIX,trialIX,:)),[-round(ssTraces{cellIX}.delay(odorIX)*1e4)+10e4 0]);
                if ~isnan(trace00(1)) && trace00(1) && any(str2num(ssTraces{cellIX}.stableResponses{2})==ssTraces{cellIX}.trialsVC00(odorIX,trialIX))
                    trace00 = conv(trace00,smoothing,'same');
                    trace00 = trace00(1:10:end);
                    trace00 = (trace00 - nanmean(trace00(1:0.8e4)));% /nanstd(trace00(1:0.8e4));
                    trace00all(counter00) = mean(trace00(integWindow));
                    trace00odor(counter00) = odorIX;
                    trace00cell(counter00) = cellIX;
                    trace00trial(counter00) = ssTraces{cellIX}.trialsVC00(odorIX,trialIX);
                    trace00odorid{counter00} = sTraces{cellIX}.odorsVCII{odorIX}(1:3);
                    counter00 = counter00 + 1;
                end
            end
            for trialIX = 1:size(sTraces{cellIX}.onsetVC70,2)
                trace70 = circshift(squeeze(sTraces{cellIX}.onsetVC70(odorIX,trialIX,:)),[-round(ssTraces{cellIX}.delay(odorIX)*1e4)+10e4 0]);
                if ~isnan(trace70(1)) && trace70(1) && any(str2num(ssTraces{cellIX}.stableResponses{1})==ssTraces{cellIX}.trialsVC70(odorIX,trialIX))
                    trace70 = conv(trace70,smoothing,'same');
                    trace70 = trace70(1:10:end);
                    trace70 = (trace70 - nanmean(trace70(1:0.8e4)));%/nanstd(trace70(1:0.8e4));
                    trace70all(counter70) = mean(trace70(integWindow));
                    trace70odor(counter70) = odorIX;
                    trace70cell(counter70) = cellIX;
                    trace70trial(counter70) = ssTraces{cellIX}.trialsVC70(odorIX,trialIX);
                    trace70odorid{counter70} = sTraces{cellIX}.odorsVCII{odorIX}(1:3);
                    counter70 = counter70 + 1;
                end
            end
        end
    end
end


rng(46e2)
% odor-specific mean currents subtract explain variance ...
trace70allExpl3 = trace70all;
trace00allExpl3 = trace00all;
trace70allExpl3x = trace70all;
trace00allExpl3x = trace00all;

odorNames = unique(trace70odorid);
odorNames2 = unique(trace00odorid);
clear odorResponseMean70 odorResponseMean00 odorResponseMean00x odorResponseMean70x
for k = 1:numel(odorNames)
    odorix = find(strcmp(trace70odorid,odorNames(k)) & trace70all < -5);
    odorResponseMean70(k) = mean(trace70all(odorix));
    odorResponseMean70x(k) = mean(trace70all(randperm(numel(trace70all),numel(odorix))));
    odorix2 = find(strcmp(trace00odorid,odorNames(k)) & trace00all > 5);
    odorResponseMean00(k) = mean(trace00all(odorix2));
    odorResponseMean00x(k) = mean(trace00all(randperm(numel(trace00all),numel(odorix))));
end
for j = 1:numel(trace70allExpl3)
    odorixx = find(strcmp(trace70odorid(j),odorNames));
    trace70allExpl3(j) = trace70allExpl3(j) - odorResponseMean70(odorixx);
    trace70allExpl3x(j) = trace70allExpl3x(j) - odorResponseMean70x(odorixx);
end
for j = 1:numel(trace00allExpl3)
    odorixx = find(strcmp(trace00odorid(j),odorNames));
    trace00allExpl3(j) = trace00allExpl3(j) - odorResponseMean00(odorixx);
    trace00allExpl3x(j) = trace00allExpl3x(j) - odorResponseMean00x(odorixx);
end



ft = fittype( 'a*x+b', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( ft );
opts.Display = 'Off';
trace70allExpl1 = trace70all;
trace00allExpl1 = trace00all;
trace70allExpl1x = trace70all;
trace00allExpl1x = trace00all;
trace70allExpl2 = trace70all;
trace00allExpl2 = trace00all;
trace70allExpl2x = trace70all;
trace00allExpl2x = trace00all;
trace70odorx = datasample(trace70odor,numel(trace70odor));
trace00odorx = datasample(trace70odor,numel(trace00odor));
warning off
clear var70explained var00explained var00 var70 var70explainedx var00explainedx var70explained2 var70explained2x var00explained2 var00explained2x var70explained3 var00explained3
for k = 1:numel(ssTraces)
    cells = find(trace70cell == k);
    if ~isempty(cells) && numel(unique(trace70odor(cells)))==3 && median(trace70all(cells)) < -0
        % tuning explains variance ..
        for jj = 1:3
            odors = cells(find(trace70odor(cells) == jj));
            odorsScrambled = cells(find(trace70odorx(cells) == jj));
            trace70allExpl1(odors) = trace70all(odors) - mean(trace70all(odors));
            trace70allExpl1x(odorsScrambled) = trace70all(odorsScrambled) - mean(trace70all(odorsScrambled));
        end
        % adaptation explains variance ..
        [xData, yData] = prepareCurveData(trace70trial(cells), trace70all(cells) );
        [fitresult, gof] = fit( xData, yData, ft, opts );
        trace70allExpl2(cells) = trace70all(cells) - fitresult(trace70trial(cells))';
        [xData, yData] = prepareCurveData(trace70trial(cells(randperm(numel(cells)))), trace70all(cells) ); % scrambled data
        [fitresult, gof] = fit( xData, yData, ft, opts );
        trace70allExpl2x(cells) = trace70all(cells) - fitresult(xData)';
        % calculate variances for this cell
        var70(k) = var(trace70all(cells));
        var70explained(k) = var(trace70allExpl1(cells));
        var70explainedx(k) = var(trace70allExpl1x(cells));
        var70explained2(k) = var(trace70allExpl2(cells));
        var70explained2x(k) = var(trace70allExpl2x(cells));
        var70explained3(k) = var(trace70allExpl3(cells));
        var70explained3x(k) = var(trace70allExpl3x(cells));
    else
        var70(k) = NaN;
        var70explained(k) = NaN;
        var70explainedx(k) = NaN;
        var70explained2(k) = NaN;
        var70explained2x(k) = NaN;
        var70explained3(k) = NaN;
        var70explained3x(k) = NaN;
    end
    cells = find(trace00cell == k);
    if ~isempty(cells) && numel(unique(trace00odor(cells)))==3 && median(trace00all(cells)) >0
        % tuning explains variance ...
        for jj = 1:3
            odors = cells(find(trace00odor(cells) == jj));
            odorsScrambled = cells(find(trace00odorx(cells) == jj));
            trace00allExpl1(odors) = trace00all(odors) - mean(trace00all(odors));
            trace00allExpl1x(odorsScrambled) = trace00all(odorsScrambled) - mean(trace00all(odorsScrambled));
        end
        % adaptation explains variance ..
        [xData, yData] = prepareCurveData(trace00trial(cells), trace00all(cells) );
        [fitresult, gof] = fit( xData, yData, ft, opts );
        trace00allExpl2(cells) = trace00all(cells) - fitresult(trace00trial(cells))';
        [xData, yData] = prepareCurveData(trace00trial(cells(randperm(numel(cells)))), trace00all(cells) ); % scrambled data
        [fitresult, gof] = fit( xData, yData, ft, opts );
        trace00allExpl2x(cells) = trace00all(cells) - fitresult(xData)';
        % calculate variances for this cell
        var00(k) = var(trace00all(cells));
        var00explained(k) = var(trace00allExpl1(cells));
        var00explainedx(k) = var(trace00allExpl1x(cells));
        var00explained2(k) = var(trace00allExpl2(cells));
        var00explained2x(k) = var(trace00allExpl2x(cells));
        var00explained3(k) = var(trace00allExpl3(cells));
        var00explained3x(k) = var(trace00allExpl3x(cells));
    else
        var00(k) = NaN;
        var00explained(k) = NaN;
        var00explainedx(k) = NaN;
        var00explained2(k) = NaN;
        var00explained2x(k) = NaN;
        var00explained3(k) = NaN;
        var00explained3x(k) = NaN;
   end
end
warning on

LL = [];
LLix = [];
% inh tuning
ixL = find(var00~=0);
LL = [LL,((var00(ixL)- var00explained(ixL))./var00(ixL))'];
LLix = [LLix,1*ones(size(var00(ixL)))];
LL = [LL,((var00(ixL)- var00explainedx(ixL))./var00(ixL))'];
LLix = [LLix,2*ones(size(var00(ixL)))];
% exc tuning
ixL = find(var70~=0);
LL = [LL,((var70(ixL)- var70explained(ixL))./var70(ixL))'];
LLix = [LLix,3*ones(size(var70(ixL)))];
LL = [LL,((var70(ixL)- var70explainedx(ixL))./var70(ixL))'];
LLix = [LLix,4*ones(size(var70(ixL)))];
% inh adaptation
ixL = find(var00~=0);
LL = [LL,((var00(ixL)- var00explained2(ixL))./var00(ixL))'];
LLix = [LLix,5*ones(size(var00(ixL)))];
LL = [LL,((var00(ixL)- var00explained2x(ixL))./var00(ixL))'];
LLix = [LLix,6*ones(size(var00(ixL)))];
% exc adaptation
ixL = find(var70~=0);
LL = [LL,((var70(ixL)- var70explained2(ixL))./var70(ixL))'];
LLix = [LLix,7*ones(size(var70(ixL)))];
LL = [LL,((var70(ixL)- var70explained2x(ixL))./var70(ixL))'];
LLix = [LLix,8*ones(size(var70(ixL)))];
% inh odor strength
ixL = find(var00~=0);
LL = [LL,((var00(ixL)- var00explained3(ixL))./var00(ixL))'];
LLix = [LLix,9*ones(size(var00(ixL)))];
% LL = [LL,((var00(ixL)- var00explained3x(ixL))./var00(ixL))'];
% LLix = [LLix,10*ones(size(var00(ixL)))];
% exc odor strength
ixL = find(var70~=0);
LL = [LL,((var70(ixL)- var70explained3(ixL))./var70(ixL))'];
LLix = [LLix,11*ones(size(var70(ixL)))];
% LL = [LL,((var70(ixL)- var70explained3x(ixL))./var70(ixL))'];
% LLix = [LLix,12*ones(size(var70(ixL)))];


figure(96), boxplot(100*(LL(:)),LLix(:),'colorgroup',[1 2 3 4 5 6 7 8 9 10 ],'notch','off','outliersize',1,'labels',{'Tuning inh','scrambled','Tuning exc','scrambled','Adaptation inh','scrambled','Adaptation exc','scrambled','Odor inh','Odor exc'});





%% point neuron model investigations (pDp; aDp, comparison with CC data)



criterion = Anatomicallocationsoverview(:,7) == 1;% | Anatomicallocationsoverview(:,3) == 1; % only Dp
% criterion = any(Anatomicallocationsoverview(:,2:8)' == 1); % any cell

integWindow = 1e4:1.15e4;
counter70 = 1; counter00 = 1; counterOdor = 1;
trace00all = zeros(numel00,1); trace70all = zeros(numel70,1); % all trials, odors, neurons
smoothing = fspecial('gaussian',[10 1],5);


Rin00 = []; base00 = []; Rin70 = []; base70 = [];
for cellIX = 114:133% 1:numel(sTraces)
    cellIX
    [onsetVC70, onsetVC00,odorsVC,dateID] = extractVC(cellIX,[],0);
    
    clear Rin base Rin_short base_short
    for k = 1:size(onsetVC70,1)
        for j = 1:size(onsetVC70,2)
            trace70 = squeeze(onsetVC70(k,j,:));
            excerpt = trace70(350001:370000);
            excerpt = reshape(excerpt,1000,20);
            baseline = excerpt(650:1000,:);
            inline = excerpt(250:490,:);
            Racc1(k,j) = 5000/(median(max(excerpt)) - median(baseline(:)));
            Racc2(k,j) = 5000/(median(baseline(:)) - median(min(excerpt)));
            Rin(k,j,:) = 5./( median(baseline) - median(inline) );
            base(k,j,:) = median(baseline) - nanmedian(trace70(2e5:end));
            Rin_short(k,j) = 5./( median(baseline(:)) - median(inline(:)) );
            base_short(k,j) = median(baseline(:)) - nanmedian(trace70(2e5:end));
        end
    end
    base70 = [base70;-base(:)];
    Rin70 = [Rin70;Rin(:)];
%     figure(29), plot(-base(:),Rin(:),'.'); hold on; %plot(-base_short(:),Rin_short(:),'r.','Markersize',18);
%     axis([-7.5 40 0 4])
    
    clear Rin base Rin_short base_short
    for k = 1:size(onsetVC00,1)
        for j = 1:size(onsetVC00,2)
            trace00 = squeeze(onsetVC00(k,j,:));
            excerpt = trace00(350001:370000);
            excerpt = reshape(excerpt,1000,20);
            baseline = excerpt(650:1000,:);
            inline = excerpt(250:490,:);
            Racc1(k,j) = 5/(median(max(excerpt)) - median(baseline(:)));
            Racc2(k,j) = 5/(median(baseline(:)) - median(min(excerpt)));
            Rin(k,j,:) = 5./( median(baseline) - median(inline) );
%             base(k,j,:) = median(baseline) - quantile(conv(trace00(2e5:end),fspecial('gaussian',[80 1],40),'same'),0.2);%nanmedian(trace00(2e5:end));
            base(k,j,:) = median(inline) - quantile(conv(trace00(2e5:end),fspecial('gaussian',[80 1],40),'same'),0.2);%nanmedian(trace00(2e5:end));
            Rin_short(k,j) = 5./( median(baseline(:)) - median(inline(:)) );
            base_short(k,j) = median(baseline(:)) - nanmedian(trace00(2e5:end));
            
        end
    end
    base00 = [base00;base(:)];
    Rin00 = [Rin00;Rin(:)];
%     figure(28), plot(base(:),Rin(:),'.'); hold on; plot(base_short(:),Rin_short(:),'r.','Markersize',18); 
%     axis([-7.5 40 0 4])

end
    
V = [Rin00,base00];
W = [Rin70,base70];
Q = [V;W];

X = Q;

X(X(:,1)>8 | X(:,1)<0) = NaN;
X(X(:,2)>100 | X(:,2)<-8) = NaN;

X1 = X(:,1);
X2 = X(:,2);


[N,C] = hist3(X,[200 200]);

figure, imagesc(C{2},C{1},conv2(N',fspecial('gaussian',[5 5],5),'same')')
axis xy


figure(28), plot(X2(:),X1(:),'.');
%     axis([-7.5 40 0 4])
    
    
    for odorIX = 1:size(sTraces{cellIX}.onsetVC70,1)
        odorCounter00 = 0;
        for trialIX = 1:size(sTraces{cellIX}.onsetVC00,2)
            trace00 = circshift(squeeze(sTraces{cellIX}.onsetVC00(odorIX,trialIX,:)),0*[-round(ssTraces{cellIX}.delay(odorIX)*1e4)+10e4 0]);
            
            excerpt = trace00(350001:370000);
            excerpt = reshape(excerpt,1000,20);
            figure, imagesc(excerpt)
            
            baseline = excerpt(650:1000,:);
            inline = excerpt(250:490,:);
            Racc1 = (median(max(excerpt)) - median(baseline(:)))/5;
            Racc2 = (median(baseline(:)) - median(min(excerpt)))/5;
            Rin = ( median(baseline(:)) - median(inline(:)) )/5;

        end
        for trialIX = 1:size(sTraces{cellIX}.onsetVC70,2)
            trace70 = circshift(squeeze(sTraces{cellIX}.onsetVC70(odorIX,trialIX,:)),[-round(ssTraces{cellIX}.delay(odorIX)*1e4)+10e4 0]);



        end
    end
end






1. determine input resistance, firing threshold and access resistance for aDp recordings

%% typical local stacks for aDp, pDp and NT >> see data

