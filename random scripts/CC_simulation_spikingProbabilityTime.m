
IXX = [97 22 22 38 38 39 39 39 70 62 70 70 86 87 97 98];
JX = [2 1 2 1 2 1 2 3 1 1 2 3 1 2 1 1];
% IXX = [67];
% JX = [1];
good_neurons = find(sum(indizes(:,[1 3]),2)==2);



load('C:\Data\rupppete\PhD\electrophysiology2016\SingleCells\OnsetDays.mat')
for IX = 1:98
    day{IX} = X_all{IX}.dateID;
end
days = unique(day);
%


nb_neurons = numel(IXX);

Trajet = zeros(nb_neurons,2e4,100);
R0_save = zeros(nb_neurons,100);
Trajet2 = zeros(nb_neurons,2e4,100);
clear SpikingMs SpikingMs2 SuperTimeStd SuperTime
for p = 1:numel(IXX)
    p
    odorIX = JX(p);
    IX = IXX(p);

    odorIndex = find(~cellfun(@isempty,strfind(datasetSingleCells{IX}.odors,X_all{IX}.odorsVC{odorIX})));
    dayIndex = find(~cellfun(@isempty,strfind(days,X_all{IX}.dateID)));
    delay = onsetZ(dayIndex,odorIndex);

    x1 = squeeze(nanmean(X_all{IX}.onsetVC70(odorIX,:,:),2))*1e-12/70e-3;
    x2 = squeeze(nanmean(X_all{IX}.onsetVC00(odorIX,:,:),2))*1e-12/70e-3;

    x1 = circshift(x1,[1e5-round(delay) 1]);
    x2 = circshift(x2,[1e5-round(delay) 1]);
%     x1 = circshift(x1,[40*uu 1]);

    x1 = x1 - median(x1(1:0.9e5));
    x2 = x2 - median(x2(1:0.9e5));
%     Cap = 20e-12;
%     R0 = 1.5e9;
    VI = -70e-3;
    VE = 0e-3;
    parallelism = 100;
    V0 = -70e-3*ones(parallelism,1);
    V = -70e-3*ones(parallelism,1);
    V2 = -70e-3*ones(parallelism,1);
    R0 = max(0.3,normrnd(2e9,0.5e9,parallelism,1)); % center, std
    Cap = 20e-12; % center, std
    dtC = 1e-4./Cap;

    VV = zeros(200000,parallelism);
    VV2 = zeros(200000,parallelism);
    for t = 1:200000
        V = V + dtC.*(  1./R0.*(V0 - V) + x2(t)*(VI - V) - x1(t)*(VE - V));
        V2 = V2 + dtC.*(  1./R0.*(V0 - V2) + 0*x2(t)*(VI - V2) - x1(t)*(VE - V2));
        VV(t,:) = V;
        VV2(t,:) = V2;
    end
    timeX = (1:2e5)/1e4;
    
    Trajet(p,:,:) = VV(1:10:end,:);
    R0_save(p,:) = R0;
    Trajet2(p,:,:) = VV2(1:10:end,:);
end
figure, imagesc(squeeze(mean(Trajet,3)))


V_thresh_center = -38; %

y = 1./(1+exp(-(Trajet-V_thresh_center/1000)/0.001)); % firing threshold around 38, width ca. 5 mV
y2 = 1./(1+exp(-(Trajet2-V_thresh_center/1000)/0.001)); % firing threshold around 38, width ca. 5 mV
figure, imagesc(squeeze(mean(y2,3)))


for hh = 1:nb_neurons
    SpikingTime(hh) = sum(sum(y(hh,10001:12000,:)))/parallelism
    SpikingTime2(hh) = sum(sum(y2(hh,10001:12000,:)))/parallelism
end

figure(4312); 
 plot([0.6+zeros(size(SpikingTime));1.4+zeros(size(SpikingTime2))],[SpikingTime;SpikingTime2],'Color',[0.9 0.9 0.9]); hold on; 
 plot(0.6+zeros(size(SpikingTime)),SpikingTime,'ko')
plot(1.4+zeros(size(SpikingTime2)),SpikingTime2,'ro'); 
plot(1.4,mean(SpikingTime2),'r.','MarkerSize',24)
plot(0.6,mean(SpikingTime),'k.','MarkerSize',24)
set(gca, 'xtick', [0.6 1.4], 'xticklabel', {'w/ inhibition', 'w/o inhibition'}, 'xlim', [0.2 1.8]);
xlim([0.2 1.8]); ylabel('Spiking probability time, first 2 sec after onset')
hold off



    V_thresh_center = -38; %

    y = 1./(1+exp(-(VV-V_thresh_center/1000)/0.001)); % firing threshold around 38, width ca. 5 mV
%         y = mean(y,2);
    SpikingMs(Rx,:) = sum(y,1)/10; % ms of probability of spiking
    y2 = 1./(1+exp(-(VV2-V_thresh_center/1000)/0.001)); % firing threshold around 38, width ca. 5 mV
%         y2 = mean(y2,2);
    SpikingMs2(Rx,:) = sum(y2,1)/10; % ms of probability of spiking
    SuperTime(counter,:) = [nanmean(SpikingMs(:)),nanmean(SpikingMs2(:))];
    SuperTimeStd(counter,:) = [nanstd(SpikingMs(:)),nanstd(SpikingMs2(:))];
    counter = counter + 1;
end
% hold off;
