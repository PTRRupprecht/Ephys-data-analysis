
IXX = [97 22 22 38 38 39 39 39 70 62 70 70 86 87 97 98];
JX = [2 1 2 1 2 1 2 3 1 1 2 3 1 2 1 1];
IXX = [67];
JX = [1];

%
load('C:\Data\rupppete\PhD\electrophysiology2016\SingleCells\OnsetDays.mat')
for IX = 1:98
    day{IX} = X_all{IX}.dateID;
end
days = unique(day);
%
% showRawData(IX)
for p = 1:numel(IXX)
    odorIX = JX(p);
    IX = IXX(p)
    
    odorIndex = find(~cellfun(@isempty,strfind(datasetSingleCells{IX}.odors,X_all{IX}.odorsVC{odorIX})));
    dayIndex = find(~cellfun(@isempty,strfind(days,X_all{IX}.dateID)));
    delay = onsetZ(dayIndex,odorIndex);

    x1 = squeeze(nanmean(X_all{IX}.onsetVC70(odorIX,:,:),2))*1e-12/70e-3;
    x2 = squeeze(nanmean(X_all{IX}.onsetVC00(odorIX,:,:),2))*1e-12/70e-3;
    
    x1 = circshift(x1,[1e5-round(delay) 1]);
    x2 = circshift(x2,[1e5-round(delay) 1]);
%     x1 = circshift(x1,[40*uu 1]);
    
    x1 = x1 - median(x1(1:1e5));
    x2 = x2 - median(x2(1:1e5));
%     figure(88); plot(x1); hold on; plot(x2)
    Cap = 20e-12;
    dtC = 1e-4/Cap;
    R0 = 1.5e9;
    VI = -70e-3;
    VE = 0e-3;
    V0 = -70e-3;
    V = -70e-3;
    V2 = -70e-3;
    VV = zeros(200000,1);
    VV2 = zeros(200000,1);
    for t = 1:200000
        V = V + dtC*(  1/R0*(V0 - V) + x2(t)*(VI - V) - x1(t)*(VE - V));
        V2 = V2 + dtC*(  1/R0*(V0 - V2) + 0*x2(t)*(VI - V2) - x1(t)*(VE - V2));
        VV(t) = V;
        VV2(t) = V2;
    end
    timeX = (1:2e5)/1e4;
%     figure(4111); plot(timeX,VV*1e3); hold on; plot(timeX,VV2*1e3,'r'); ylim([-80 30])
%     VVV(:,p) = VV;
%     VVV2(:,p) = VV2;
    VV(VV>-0.038) = 0;
    VV2(VV2>-0.038) = 0;
    
    drawnow;
    figure(41892), subplot(3,1,1); plot(timeX,VV*1e3); hold on; plot(timeX,VV2*1e3,'r'); hold off; ylim([-80 30])
    xlabel('time [sec]'); ylabel('membrane potential (simulation) [mV]')
    subplot(3,1,2); plot(timeX,smooth(-x1,100)); ylabel('Excitation [Sievert]'); xlabel('time [sec]')
    subplot(3,1,3); plot(timeX,smooth(x2,100)); ylabel('Inhibition [Sievert]'); xlabel('time [sec]')
    pause(1)
end
% hold off;
