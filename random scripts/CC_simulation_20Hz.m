
IXX = [97 22 22 38 38 39 39 39 70 62 70 70 86 87 97 98];
JX = [2 1 2 1 2 1 2 3 1 1 2 3 1 2 1 1];
IXX = [65 67 70 86 96 ];
JX = [1 1 1 2 2];

%
load('C:\Data\rupppete\PhD\electrophysiology2016\SingleCells\OnsetDays.mat')
for IX = 1:98
    day{IX} = X_all{IX}.dateID;
end
days = unique(day);
%

counter = 1;
clear SpikingMs SpikingMs2 SuperTimeStd SuperTime
superVoltage = zeros(20,90,100);
for p = 1:numel(IXX)
    p
    
    for reps = 1:20
        odorIX = JX(p);
        IX = IXX(p);
        odorIndex = find(~cellfun(@isempty,strfind(datasetSingleCells{IX}.odors,X_all{IX}.odorsVC{odorIX})));
        dayIndex = find(~cellfun(@isempty,strfind(days,X_all{IX}.dateID)));
        delay = onsetZ(dayIndex,odorIndex);
        odorchoice1 = randi(size(X_all{IX}.onsetVC70(odorIX,:,:),2));
        odorchoice2 = randi(size(X_all{IX}.onsetVC00(odorIX,:,:),2));

        x1 = squeeze(nanmean(X_all{IX}.onsetVC70(odorIX,odorchoice1,:),2))*1e-12/70e-3;
        x2 = squeeze(nanmean(X_all{IX}.onsetVC00(odorIX,odorchoice2,:),2))*1e-12/70e-3;

        x1 = circshift(x1,[1e5-round(delay) 1]);
        x2 = circshift(x2,[1e5-round(delay) 1]);
        x1 = x1 - median(x1(1:1e5));
        x2 = x2 - median(x2(1:1e5));
        clear x1dd x2dd
        for dd = 1:90
            x1dd(dd,:) = circshift(x1,[10*dd-450 1]);
            x2dd(dd,:) = x2;
        end
        x1dd = x1dd(:,70001:150000);
        x2dd = x2dd(:,70001:150000);
        
    %     Cap = 20e-12;
    %     R0 = 1.5e9;
        VI = -70e-3;
        VE = 0e-3;
        parallelism = 40;
        V0 = -70e-3*ones(parallelism,1);
        V = -70e-3*ones(parallelism,90);
        V2 = -70e-3*ones(parallelism,1);
        R0 = max(0.3,normrnd(2.5e9,0.5e9,parallelism,1)); % center, std
        Cap = 20e-12;
        dtC = 1e-4./Cap;

        VV = zeros(80000,parallelism,90);
        for t = 1:80000
            if mod(t,20000)==0; disp(t); end
            for dd = 1:90
                V(:,dd) = V(:,dd) + dtC*(  1./R0.*(V0 - V(:,dd)) + x2dd(dd,t)*(VI - V(:,dd)) - x1dd(dd,t)*(VE - V(:,dd)));
                VV(t,:,dd) = V(:,dd);
            end
        end

        threshes = 1-1./exp(linspace(log(1/(1-0.95)),log(1/(1-0.9999)),20));

    %     ff = squeeze(mean(max(VV(1e5:1.4e5,:,:),[],1),2));
        for kk = 1:numel(threshes)
%             kk
            superVoltage(kk,:,counter) = squeeze(mean(quantile(VV(2.5e4:6e4,:,:),threshes(kk),1),2));
        end
        choicesodor1(counter) = odorchoice1;
        choicesodor2(counter) = odorchoice2;
        
        counter = counter + 1;
    end
end
save('temp.mat','superVoltage','choicesodor1','choicesodor2','IXX','JX');
% hold off;


%% get corresponding traces for visualization
IXX = [65 67 70 70 86  ];
JX = [1 1 1 1 2 ];
odor1 = choicesodor1([5 25 48 54 69]);
odor2 = choicesodor2([5 25 48 54 69]);


p = 1
odorIX = JX(p);
IX = IXX(p);
odorIndex = find(~cellfun(@isempty,strfind(datasetSingleCells{IX}.odors,X_all{IX}.odorsVC{odorIX})));
dayIndex = find(~cellfun(@isempty,strfind(days,X_all{IX}.dateID)));
delay = onsetZ(dayIndex,odorIndex);
odorchoice1 = odor1(p);
odorchoice2 = odor2(p);

x1 = squeeze(nanmean(X_all{IX}.onsetVC70(odorIX,odorchoice1,:),2))*1e-12/70e-3;
x2 = squeeze(nanmean(X_all{IX}.onsetVC00(odorIX,odorchoice2,:),2))*1e-12/70e-3;

x1 = circshift(x1,[1e5-round(delay) 1]);
x2 = circshift(x2,[1e5-round(delay) 1]);
x1 = x1 - median(x1(1:1e5));
x2 = x2 - median(x2(1:1e5));
clear x1dd x2dd
for dd = 1:90
    x1dd(dd,:) = circshift(x1,[10*dd-450 0]);
    x2dd(dd,:) = x2;
end
x1dd = x1dd(:,70001:150000);
x2dd = x2dd(:,70001:150000);

%     Cap = 20e-12;
%     R0 = 1.5e9;
VI = -70e-3;
VE = 0e-3;
parallelism = 40;
V0 = -70e-3*ones(parallelism,1);
V = -70e-3*ones(parallelism,90);
V2 = -70e-3*ones(parallelism,1);
R0 = max(0.3,normrnd(2.5e9,0.5e9,parallelism,1)); % center, std
Cap = 20e-12;
dtC = 1e-4./Cap;

VV = zeros(80000,parallelism,90);
for t = 1:80000
    if mod(t,20000)==0; disp(t); end
    for dd = 1:90
        V(:,dd) = V(:,dd) + dtC*(  1./R0.*(V0 - V(:,dd)) + x2dd(dd,t)*(VI - V(:,dd)) - x1dd(dd,t)*(VE - V(:,dd)));
        VV(t,:,dd) = V(:,dd);
    end
end
figure(4), subplot(2,3,3);
plot((1:80000)/1e4-3,circshift(mean(VV(:,:,82)')*1e3,[0 -215])); hold on; plot((1:80000)/1e4-3,circshift(mean(VV(:,:,61)')*1e3,[-0 0]),'r')
xlim([-1 3]); ylim([-80 -30])



