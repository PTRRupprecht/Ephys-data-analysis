



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





counterCell = 1;
clear AUCs2 AUCs2multi AUCs2permute AUCspermute
for k = 1:numel(sTraces)
    trials = find(trace70cell == k);
    trials2 = find(trace00cell == k);
    if ~isempty(trials) &&  ~isempty(trials2)
        AUCs = trace70all(trials);
        AUCs2 = trace00all(trials2);
%         for pp = 1:100
%             AUCspermute{pp} = AUCs(randperm(numel(AUCs(:))));
%         end
        odors = trace70odor(trials);
        odors2 = trace00odor(trials2);
        if numel(unique(odors)) == 3 && numel(unique(odors2)) == 3
            for jj = 1:3
                trialsX = find(odors == jj);
                trialsX2 = find(odors2 == jj);
                AUCsX(counterCell,jj) = mean(AUCs(trialsX));
                AUCsX2(counterCell,jj) = mean(AUCs2(trialsX2));
%                 for pp = 1:100
%                     AUCs2permute{pp}(counterCell,jj) = mean(AUCspermute{pp}(trialsX));
%                 end
%                 AUCs2multi(counterCell,jj) = numel(trialsX);
            end
            counterCell = counterCell + 1;
        end
    end
end

clear IX
for j = 1:size(AUCsX,1)
    [~,IXX] = sort(AUCsX2(j,:),'descend');
    [~,IXX2] = sort(AUCsX(j,:));
%     IX(:,j) = IXX;
    
    AUCsXX(j,:) = AUCsX(j,IXX);
    AUCsXX2(j,:) = AUCsX2(j,IXX);
    AUCsXY(j,:) = AUCsX(j,IXX2);
    AUCsXY2(j,:) = AUCsX2(j,IXX2);
end


figure(32), plot(-mean(AUCsXX)','b'); hold on; plot(mean(AUCsXX2)','r');
 plot(-mean(AUCsXY),'--'); hold on; plot(mean(AUCsXY2),'--r');hold off

    [~,IX] = sort(AUCsX');

figure, imagesc(AUCsX(IX'))

