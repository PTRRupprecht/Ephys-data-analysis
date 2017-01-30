

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


for j = 1:size(Anatomicallocationsoverview,1)
    if Anatomicallocationsoverview(j,2) == 1 || Anatomicallocationsoverview(j,4) == 1
        Anatomicallocationsoverview(j,3) = NaN;
    elseif Anatomicallocationsoverview(j,7) == 1
        Anatomicallocationsoverview(j,6) = NaN;
        Anatomicallocationsoverview(j,5) = NaN;
    end
end
        
[anatomicalColorsY,anatomicalColorsX] = find(~isnan(Anatomicallocationsoverview(:,2:end)));
NeuronViewerMap(anatomicalColorsX,anatomicalColorsY,513524)





%% real data analysis

numel00 = 0; numel70 = 0;
for cellIX = 1:numel(sTraces)
    for odorIX = 1:size(sTraces{cellIX}.onsetVC70,1)
        numel00 = numel00+size(sTraces{cellIX}.onsetVC00,2);
        numel70 = numel70+size(sTraces{cellIX}.onsetVC70,2);
    end
end

criterion = Anatomicallocationsoverview(:,7) == 1;% | Anatomicallocationsoverview(:,3) == 1; % only Dp
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
% rankingResponse = nanmean(chosen70odor(:,2.1e4:2.5e4),2);
[~,XI] = sort(rankingResponse,'ascend');
figure(2),subplot(1,2,1); imagesc([0 60],[],chosen70odor(XI,:),[-2 2]); xlim([5 20]);
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
AUC_00 = mean(chosen00odor(:,integWindow)');%*numel(integWindow)/1e3;
AUC_70 = mean(chosen70odor(:,integWindow)');%*numel(integWindow)/1e3;
integWindow = 1.2e4:1.9e4;
AUC_002 = mean(chosen00odor(:,integWindow)');%*numel(integWindow)/1e3;
AUC_702 = mean(chosen70odor(:,integWindow)');%*numel(integWindow)/1e3;

figure(10); 
for k = 1:numel(AUC_00)
    plot([AUC_00(k) ],-[AUC_70(k) ],'.g','Markersize',16); hold on
    plot([ AUC_002(k)],-[AUC_702(k)],'.b','Markersize',16); hold on
    plot([AUC_00(k) AUC_002(k)],-[AUC_70(k) AUC_702(k)],'-','Markersize',16)
end
plot(mean(AUC_00),-mean(AUC_70),'.g','Markersize',36); hold on
plot(mean(AUC_002),-mean(AUC_702),'.b','Markersize',36); hold on
plot([mean(AUC_00) mean(AUC_002)],[-mean(AUC_70) -mean(AUC_702)],'-','LineWidth',3)
axis([-5 45 -5 37.5])



%% phase plot for exc. and inhibition (trajectories in EPSC/IPSC space)

cmap = jet(150);
counter = 0;
ch70All = [];%zeros(1501,1);
ch00All = [];%zeros(1501,1);
figure(54),
for k = 1:size(chosen00odor,1)
    ch00 = smooth(chosen00odor(k,1e4:2e4),250);
    ch70 = -smooth(chosen70odor(k,1e4:2e4),250);
    if max(ch00) > 5 && max(ch70) > 5
        counter = counter + 1;
%         ch00 = (ch00 - quantile(ch00(1:200),0.005))/(quantile(ch00,0.995) - quantile(ch00(1:200),0.005));
        ch00 = (ch00 - mean(ch00(1:10)))/(max(ch00) - mean(ch00(1:10)));
        ch00All = [ch00All, ch00];
        ch70 = (ch70 - mean(ch70(1:10)))/(max(ch70) - mean(ch70(1:10)));
%         ch70 = (ch70 - quantile(ch70(1:200),0.005))/(quantile(ch70,0.995) - quantile(ch70(1:200),0.005));
        ch70All = [ch70All, ch70];
        for j = 1:10:1350
            plot(smooth(ch70((j):(j+51)),10),smooth(ch00((j):(j+51)),10),'Color',cmap(min(150,ceil(j/5)),:)); hold on;
        end
    end
end;
% ch70All = median(ch70All');
% ch00All = median(ch00All');
ch70All = mean(ch70All');
ch00All = mean(ch00All');
plot(ch70All(1:1350),ch00All(1:1350),'k','LineWidth',7); hold on;
for j = 1:10:1350
    plot(ch70All((j):(j+51)),ch00All((j):(j+51)),'Color',cmap(min(150,ceil(j/5)),:),'LineWidth',5); hold on;
end
hold off
axis([-.1 1.0 -0.1 1.0])
set(gca,'FontSize',13)
% colorbar of time
figure, imagesc([zeros(24),ones(24)],[0 1250])
set(gca,'FontSize',13)


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
% integWindow = 1.2e4:2e4;
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
                    trace00odorid{counter00} = sTraces{cellIX}.odorsVC{odorIX}(1:3);
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
                    trace70odorid{counter70} = sTraces{cellIX}.odorsVC{odorIX}(1:3);
                    counter70 = counter70 + 1;
                end
            end
        end
    end
end


rng(46e2)
% odor-specific mean currents subtraction-explained variance ...
trace70allExpl3 = trace70all;
trace00allExpl3 = trace00all;
trace70allExpl3x = trace70all;
trace00allExpl3x = trace00all;

odorNames = unique(trace70odorid);
odorNames2 = unique(trace00odorid);
clear odorResponseMean70 odorResponseMean00 odorResponseMean00x odorResponseMean70x
for k = 1:numel(odorNames)
    odorix = find(strcmp(trace70odorid,odorNames(k)));% & trace70all < -5);
    odorResponseMean70(k) = mean(trace70all(odorix));
    odorResponseMean70x(k) = mean(trace70all(randperm(numel(trace70all),numel(odorix))));
    odorix2 = find(strcmp(trace00odorid,odorNames(k)));% & trace00all > 5);
    odorResponseMean00(k) = mean(trace00all(odorix2));
    odorResponseMean00x(k) = mean(trace00all(randperm(numel(trace00all),numel(odorix))));
end

figure(44); hold on;
for k = 1:numel(odorNames)
    plot(-odorResponseMean70(k),odorResponseMean00(k),'.','MarkerSize',16);
    if k == 2
        text(-odorResponseMean70(k)+0.15,odorResponseMean00(k)+0.3,odorNames(k));
    elseif k == 3
        text(-odorResponseMean70(k)+0.15,odorResponseMean00(k)-0.3,odorNames(k));
    else
        text(-odorResponseMean70(k)+0.15,odorResponseMean00(k),odorNames(k));
    end
end
axis([0 12 0 15])

for j = 1:numel(trace70allExpl3)
    odorixx = find(strcmp(trace70odorid(j),odorNames));
    trace70allExpl3(j) = trace70allExpl3(j)./odorResponseMean70(odorixx)*mean(trace70all);
    trace70allExpl3x(j) = trace70allExpl3x(j)./odorResponseMean70x(odorixx)*mean(trace70all);
end
for j = 1:numel(trace00allExpl3)
    odorixx = find(strcmp(trace00odorid(j),odorNames));
    trace00allExpl3(j) = trace00allExpl3(j)./odorResponseMean00(odorixx)*mean(trace00all);
    trace00allExpl3x(j) = trace00allExpl3x(j)./odorResponseMean00x(odorixx)*mean(trace00all);
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
% LL = [LL,((var00explained3x(ixL)- var00explained3(ixL))./var00explained3x(ixL))'];
LLix = [LLix,9*ones(size(var00(ixL)))];
LL = [LL,((var00(ixL)- var00explained3x(ixL))./var00(ixL))'];
LLix = [LLix,10*ones(size(var00(ixL)))];
% exc odor strength
ixL = find(var70~=0);
LL = [LL,((var70(ixL)- var70explained3(ixL))./var70(ixL))'];
% LL = [LL,((var70explained3x(ixL)- var70explained3(ixL))./var70explained3x(ixL))'];
LLix = [LLix,11*ones(size(var70(ixL)))];
LL = [LL,((var70(ixL)- var70explained3x(ixL))./var70(ixL))'];
LLix = [LLix,12*ones(size(var70(ixL)))];


figure(96), boxplot(100*(LL(:)),LLix(:),'colorgroup',[1 2 3 4 5 6 7 8 9 10 11 12 ],'notch','off','outliersize',1,'labels',{'Tuning inh','scrambled','Tuning exc','scrambled','Adaptation inh','scrambled','Adaptation exc','scrambled','Odor inh','Odor inhx','Odor exc','Odor excx'});




LL = [];
LLix = [];
% inh tuning
ixL = find(var00~=0);
LL = [LL,((var00(ixL)- var00explained(ixL))./var00(ixL))'];
LLix = [LLix,1*ones(size(var00(ixL)))];
% inh odor strength
ixL = find(var00~=0);
LL = [LL,((var00(ixL)- var00explained3(ixL))./var00(ixL))'];
LLix = [LLix,2*ones(size(var00(ixL)))];
% exc tuning
ixL = find(var70~=0);
LL = [LL,((var70(ixL)- var70explained(ixL))./var70(ixL))'];
LLix = [LLix,3*ones(size(var70(ixL)))];
% exc odor strength
ixL = find(var70~=0);
LL = [LL,((var70(ixL)- var70explained3(ixL))./var70(ixL))'];
LLix = [LLix,4*ones(size(var70(ixL)))];

figure(98), boxplot(100*(LL(:)),LLix(:),'colorgroup',[1 2 3 4],'notch','off','outliersize',1,'labels',{'Tuning inh','Odor inh','Tuning exc','Odor exc'},'plotstyle','compact','labelorientation','horizontal');
ylim([-95 105]); hold on;
for k = 1:size(LL,1)
    shift = +0.15*(rand()-0.5);
    plot([1.25 1.75 ]+shift,LL(k,[1 2])*100,'Color',[0.34 0.34 0.34],'LineWidth',2);
    plot([1.25 1.75 ]+shift,LL(k,[1 2])*100,'.','Color',[0.6 0.6 0.6],'MarkerSize',22);
    plot(2+[1.25 1.75 ]+shift,LL(k,[ 3 4])*100,'Color',[0.34 0.34 0.34],'LineWidth',2);
    plot(2+[1.25 1.75 ]+shift,LL(k,[ 3 4])*100,'.','Color',[0.6 0.6 0.6],'MarkerSize',22);
end
hold off;

% does inhibitory tuning correlate with excitatory tuning? - No evidence found.
LL = [];
LLix = [];
% inh tuning
ixL = find(var00~=0);
LL = [LL,((var00(ixL)- var00explained(ixL))./var00(ixL))'];
LLix = [LLix,1*ones(size(var00(ixL)))];
% exc tuning
ixL = find(var70~=0);
LL = [LL,((var70(ixL)- var70explained(ixL))./var70(ixL))'];
LLix = [LLix,2*ones(size(var70(ixL)))];

figure(99), boxplot(100*(LL(:)),LLix(:),'colorgroup',[1 2 ],'notch','off','outliersize',1,'labels',{'Tuning inh','Tuning exc'});
ylim([-5 105]); hold on;
for k = 1:size(LL,1)
    shift = +0.15*(rand()-0.5);
    plot([1.25 1.75 ]+shift,LL(k,:)*100,'Color',[0.6 0.6 0.6],'LineWidth',2);
    plot([1.25 1.75]+shift,LL(k,:)*100,'.','Color',[0.34 0.34 0.34],'MarkerSize',22);
end
hold off;

chx = find(~isnan(LL(:,2)) & ~isnan(LL(:,1)));

corr(LL(chx,1),LL(chx,2))



%% point neuron model investigations (pDp; aDp, comparison with CC data)



%% choose parameters based on data
% sigma for pDp, NT, aDp; then mu for same; both of the normal distribution
% underlying the lognormal distribution
zTaum = [-2.9900    0.2012;
        -2.7913    0.2624;
        -2.9984    0.3955]'; % pDp, NT (low stat), aDp; mu and sigma
zRm = [-0.0876 0.3697 0.2379;
    0.3373 0.2709 0.2728]; % in GOhm
zCm = [4.0337 3.7091 3.681;
    0.2114 0.3671 0.4036]; % in pF
for j = 1:3
    zCm_med(j) = median(lognrnd(zCm(1,j),zCm(2,j),[100000 1]));
end
for j = 1:3
    zTaum_med(j) = median(lognrnd(zTaum(1,j),zTaum(2,j),[100000 1]));
end
% Vthresh obeys a normal distribution
zVthresh = [ 38 2.5];

% use lognormal : lognrnd(mu,sigma)
% use normal    : normrnd(mu,sigma)

clear trace00odorid cell_traces cell_tracesDisInh cell_traces_shuffled_inh cell_tracesIE cell_tracesII

for areaIX = 1:1; % 1=pDp, 2=NT, 3=aDp

%% choose set of neurons
template = [2 4 7]; % pDp, NT, aDp
areaIXX = template(areaIX);
criterion = Anatomicallocationsoverview(:,areaIXX) == 1;


%% load data

max_rep = 10;


for cellIX = 67%1:147%numel(sTraces)
    [areaIX,cellIX]
    if criterion(cellIX)
        cell_traces{cellIX} = NaN*ones(max_rep,2e3,3);
        cell_tracesII{cellIX} = NaN*ones(max_rep,2e3,3);
        cell_tracesIE{cellIX} = NaN*ones(max_rep,2e3,3);
        cell_tracesDisInh{cellIX} = NaN*ones(max_rep,2e3,3);
        cell_traces_shuffled_inh{cellIX} = NaN*ones(max_rep,2e3,3);
        for odorIX = 1:size(sTraces{cellIX}.onsetVC00,1)
            for rep = 1:max_rep
%                 Vthresh = normrnd(zVthresh(1),zVthresh(2));
                Cm = zCm_med(areaIX);
                Rm = max(0.3,1e3*lognrnd(zTaum(1,areaIX),zTaum(2,areaIX))/Cm);
                
                continueX = 0; counterX = 0;
                while continueX == 0
                    n00 = randi(size(sTraces{cellIX}.onsetVC00,2));
                    n70 = randi(size(sTraces{cellIX}.onsetVC70,2));
                    trace00 = circshift(squeeze(sTraces{cellIX}.onsetVC00(odorIX,n00,:)),[-round(ssTraces{cellIX}.delay(odorIX)*1e4)+10e4 0]);
                    odorIX_shuffled = randi(size(sTraces{cellIX}.onsetVC00,1));
                    trace00_shuffled = circshift(squeeze(sTraces{cellIX}.onsetVC00(odorIX_shuffled,n00,:)),[-round(ssTraces{cellIX}.delay(odorIX_shuffled)*1e4)+10e4 0]);
                    trace70 = circshift(squeeze(sTraces{cellIX}.onsetVC70(odorIX,n70,:)),[-round(ssTraces{cellIX}.delay(odorIX)*1e4)+10e4 0]);
                    trace00 = trace00(1:200000);
                    trace70 = trace70(1:200000);
                    if ~isnan(trace00(1)) && trace00(1) && ~isnan(trace70(1)) && trace70(1) && any(str2num(ssTraces{cellIX}.stableResponses{2})==ssTraces{cellIX}.trialsVC00(odorIX,n00)) && any(str2num(ssTraces{cellIX}.stableResponses{1})==ssTraces{cellIX}.trialsVC70(odorIX,n70))
                        
                        trace00 = (trace00 - quantile(trace00(1:0.8e5),0.33))/70e-3*1e-12;
                        trace00_shuffled = (trace00_shuffled - quantile(trace00_shuffled(1:0.8e5),0.33))/70e-3*1e-12;
                        trace70 = (trace70 - quantile(trace70(1:0.8e5),0.66))/70e-3*1e-12;
                        
                        VI = -70e-3;
                        VE = 0e-3;
                        V0 = -70e-3;
                        V = -70e-3;
                        V2 = -70e-3;
                        V3 = -70e-3;
                        R0 = 1e9*Rm; % center, std
                        Cap = Cm*1e-12; % center, std
                        dtC = 1e-4./Cap;

                        VV  = zeros(200000,1);
                        IE  = zeros(200000,1);
                        II  = zeros(200000,1);
                        VV2 = zeros(200000,1);
                        VV3 = zeros(200000,1);
                        for t = 1:200000
                            V = V + dtC.*(  1./R0.*(V0 - V) + trace00(t)*(VI - V) - trace70(t)*(VE - V));
                            IE(t) = - trace70(t)*(VE - V);
                            II(t) = + trace00(t)*(VI - V);
                            V2 = V2 + dtC.*(  1./R0.*(V0 - V2) + 0*trace00(t)*(VI - V2) - trace70(t)*(VE - V2));
                            V3 = V3 + dtC.*(  1./R0.*(V0 - V3) + trace00_shuffled(t)*(VI - V3) - trace70(t)*(VE - V3));
                            VV(t) = V;
                            VV2(t) = V2;
                            VV3(t) = V3;
                        end
                        cell_traces{cellIX}(rep,:,odorIX) = VV(1:100:end);
                        cell_traces_shuffled_inh{cellIX}(rep,:,odorIX) = VV3(1:100:end);
                        cell_tracesIE{cellIX}(rep,:,odorIX) = IE(1:100:end);
                        cell_tracesII{cellIX}(rep,:,odorIX) = II(1:100:end);
                        cell_tracesDisInh{cellIX}(rep,:,odorIX) = VV2(1:100:end);
                        continueX = 1;
                    else
                        counterX = counterX + 1;
                    end
                    if counterX > 25 % give up on this repetition
                         continueX = 1;
                    end
                end
                
            end
            trace00odorid{cellIX,odorIX} = sTraces{cellIX}.odorsVC{odorIX}(1:3);
        end
    end
end

end

cellIX = 86;
figure(3), plot(cell_tracesIE{cellIX}(:,1:1:end,2)')
hold on, plot(cell_tracesII{cellIX}(:,1:1:end,2)')
figure(5), plot(cell_tracesDisInh{cellIX}(:,1:1:end,2)','b')
hold on, plot(cell_traces{cellIX}(:,1:1:end,2)','k')


% avg over MAT, avg over FiringProb
% Probabilities for each cell-odor-pair w/ and w/o




k = 1;
clear FiringX FiringZ FiringXstd FiringZstd FiringProbW FiringProbWO
criterion = Anatomicallocationsoverview(:,2) == 1;
for cellIX = 1:147
    if criterion(cellIX) == 1
        for j = 1:3
            MAT = cell_traces{cellIX}(:,1:1:end,j);
            FiringProb = MAT;%1./(1+exp(-(MAT+42/1000)/0.001)); % firing threshold around 38, width ca. 5 mV
            FiringProbW(k,:) = mean(FiringProb);
            FiringX(k) = mean(mean(FiringProb(:,1001:1100)));
            FiringXstd(k) = std(mean(FiringProb(:,1001:1100)));
            
            MAT = cell_tracesDisInh{cellIX}(:,1:1:end,j);
            FiringProb =  MAT;%1./(1+exp(-(MAT+42/1000)/0.001)); % firing threshold around 38, width ca. 5 mV
            FiringProbWO(k,:) = mean(FiringProb);
            FiringZ(k) = mean(mean(FiringProb(:,1001:1100)));
            FiringZstd(k) = std(mean(FiringProb(:,1001:1100)));
            k = k + 1;
        end
    end
end



figure(41);
subplot(1,2,1); hold on;
for k = 1:numel(FiringZ)
    offset = rand()*0.10;
    plot(1+offset,FiringZ(k)*1000,'k.','MarkerSize',22)
    plot(2+offset,FiringX(k)*1000,'k.','MarkerSize',22)
    plot([1 2]+offset,[FiringZ(k) FiringX(k)]*1000,'k')
end
hold off; set(gca,'XTick',[1 2 ]); set(gca,'XTickLabel',{'without' 'with'});
ylabel('Average time above threshold during onset response [ms]')
xlim([0.8 2.2])

subplot(1,2,2); hold on;

plot((1:2000)/100-10,nanmean(FiringProbWO,1)*1000,'r'); hold on
plot((1:2000)/100-10,nanmean(FiringProbW,1)*1000,'k'); hold off;
ylabel('Spiking probability, averaged across all cells')
xlabel('time [sec]'); xlim([-0.5 2])

% FiringProbWX = FiringProbW;
% FiringProbWOX = FiringProbWO;
% FiringProbWX([83 91 92],:) = NaN;
% FiringProbWOX([83 91 92],:) = NaN;
% subplot(1,3,3);
% plot((1:2000)/100-10,nanmean(FiringProbWOX,1),'r'); hold on
% plot((1:2000)/100-10,nanmean(FiringProbWX,1),'k'); hold off;
% ylabel('Spiking probability, averaged across all cells')
% xlabel('time [sec]'); xlim([-0.5 2])



%% tuning a la Thomas for voltages (plus: scrambled for control)
% select region ('criterion'), inh/exc ('traceXXall' etc.) and the
% corresponding plotting modality ('order' etc., and subplot panel) and
% possiblythe integration window of interest ('integWindow')


criterion = Anatomicallocationsoverview(:,2) == 1;% | Anatomicallocationsoverview(:,3) == 1; % only Dp
% criterion = any(Anatomicallocationsoverview(:,2:8)' == 1); % any cell

clear MeanVoltage_shuffledInhibition MeanVoltage
counter = 1;
for k = 1:numel(cell_traces)
    if ~isempty(trace00odorid{k,3})
        MeanVoltage(counter,:) = squeeze(nanmean(nanmean(cell_traces{k}(:,1001:1100,:),1),2));
        [MeanVoltage(counter,:),XXI] = sort(MeanVoltage(counter,:),'descend');
        MeanVoltage_shuffledInhibition(counter,:) = squeeze(nanmean(nanmean(cell_traces_shuffled_inh{k}(:,1001:1100,:),1),2));
        MeanVoltage_shuffledInhibition(counter,:) = MeanVoltage_shuffledInhibition(counter,XXI);
        counter = counter + 1;
    end
end

range(nanmean(MeanVoltage-MeanVoltage_shuffledInhibition))/(max(nanmean(MeanVoltage_shuffledInhibition)) - min(nanmean(MeanVoltage_shuffledInhibition)))


figure(23), hold on;
plot((-70:-35)+70,(-70:-35)+70,'k');
plot(MeanVoltage_shuffledInhibition*1000+70,MeanVoltage*1000+70,'.');
hold off; axis([0 20 0 20])

figure(24);% plot((MeanVoltage-MeanVoltage_shuffledInhibition)'*1000,'-')
hold on;
bar(nanmean(MeanVoltage-MeanVoltage_shuffledInhibition)'*1000,'k','LineWidth',5)
hold off;


figure(25); bar([nanmean(MeanVoltage);nanmean(MeanVoltage_shuffledInhibition)]'*1000+70);

hold on; plot([nanmean(MeanVoltage);nanmean(MeanVoltage_shuffledInhibition)]'*1000+70);



