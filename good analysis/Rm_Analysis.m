
AAAB : [5 6 7],[8 9 10 14 15 16],[17 18 19]

AAAC :  [5 6 7 8],[9 10 11 14 15 16 17]

AAAD : [4 5 6 7],[8 9 (9)]


j = 3
% inhTrials = [5 6 7];
excTrials = [8 9 10 14 15 16];
excTrials = [9 10 11 14 15 16 17];
excTrials = [8 9 ];
excTrials = [17 18 19]

FileList = dir('**.xsg');
clear INH EXC
counter = 1;
for i = inhTrials
    load(FileList(i).name,'-mat');
    A = data.ephys.trace_1;
    INH(:,counter) = A;
    counter = counter + 1;
end
counter = 1;
for i = excTrials
    load(FileList(i).name,'-mat');
    A = data.ephys.trace_1;
    EXC(:,counter) = A;
    counter = counter + 1;
end


excMean = mean(EXC');
 
excMean = reshape(excMean,1000,380);

pixelz = 1;
baseline = (conv2(excMean(600:990,2:end),fspecial('gaussian',[1 pixelz],3),'same') + conv2(excMean(600:990,1:end-1),fspecial('gaussian',[1 pixelz],3),'same'))/2;
inline = conv2(excMean(100:490,2:end),fspecial('gaussian',[1 pixelz],3),'same');


Rin = 5./( median(baseline) - median(inline) );


IX = 150

tttrace = median(excMean');

Rs1 = 5e-3./(median(tttrace)-min(tttrace))*1e12/1e6;
Rs2 = - 5e-3./(median(tttrace)-max(tttrace))*1e12/1e6;

biophysics{IX}.Rin = Rin*1e9;
biophysics{IX}.Rs = [Rs1,Rs2];

Y = median(excMean');
X = (0:numel(Y)-1)/10;
[xData, yData] = prepareCurveData( X(2:500), Y(2:500) );
ft = fittype( 'a*exp(-x/b)+d*exp(-x/e)+c', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( ft ); opts.Display = 'Off';
opts.Lower = [-Inf -Inf -Inf -Inf -Inf];
opts.StartPoint = [0.107324056896322 20 0.374643668266154 -22 1];
opts.Upper = [Inf Inf Inf Inf Inf];
[fitresult, gof] = fit( xData, yData, ft, opts );

tauVC = max(fitresult.b,fitresult.e)

tauVC/mean([Rs1 Rs2])*1e9

biophysics{IX}.Cm = tauVC*1e-3/mean([Rs1 Rs2]);
biophysics{IX}.tau = tauVC*1e-3/mean([Rs1 Rs2])*mean(Rin(:))*1e9;


% figure(3), plot(mean(Rin),'.'); 




%% treat neurons automatically that have test pulses (sec 35 start)

Rin00 = []; base00 = []; Rin70 = []; base70 = [];
for cellIX = 115:147

    [onsetVC70, onsetVC00,odorsVC,dateID] = extractVC(cellIX,[],0);

    clear Rin base Rin_short base_short Racc1 Racc2 tauVC
    for k = 1:size(onsetVC70,1)
        for j = 1:size(onsetVC70,2)
            trace70 = squeeze(onsetVC70(k,j,:));
            excerpt = trace70(350001:370000);
            excerpt = reshape(excerpt,1000,20);
            baseline = excerpt(650:1000,:);
            inline = excerpt(250:490,:);
            Racc1(k,j) = 5000/(median(max(excerpt)) - median(baseline(:)));
            Racc2(k,j) = 5000/(median(baseline(:)) - median(min(excerpt)));

            Y = median(excerpt');
            X = (0:numel(Y)-1)/10;
            diffX = Y - circshift(Y,[0 -1]); startIDX = find(diffX<0,1,'first');
            [xData, yData] = prepareCurveData( X(startIDX:500), Y(startIDX:500) );
            ft = fittype( 'a*exp(-x/b)+d*exp(-x/e)+c', 'independent', 'x', 'dependent', 'y' );
            opts = fitoptions( ft ); opts.Display = 'Off';
            opts.Lower = [-Inf -Inf -Inf -Inf -Inf];
            opts.StartPoint = [-30 5 -22 -40 0.1];
            opts.Upper = [Inf Inf Inf Inf Inf];
            try
                [fitresult, gof,outputX] = fit( xData, yData, ft, opts );

                tauVC(j,k) = max(fitresult.b,fitresult.e)
            catch
                tauVC(j,k) = NaN;
            end

    %         Rin(k,j,:) = 5./( median(baseline) - median(inline) );
    %         base(k,j,:) = median(baseline) - nanmedian(trace70(2e5:end));
            Rin_short(k,j) = 5./( median(baseline(:)) - median(inline(:)) );
        end
    end
    
    biophysics{cellIX}.Rin = Rin_short*1e9;
    biophysics{cellIX}.Rs = [Racc1;Racc2];
    biophysics{cellIX}.Cm = tauVC.*1e-3./mean(cat(3,Racc1,Racc2),3)'/1e6;
    biophysics{cellIX}.tau = tauVC.*1e-3./mean(cat(3,Racc1,Racc2),3)'.*Rin_short'*1e9/1e6;



end

% post-analysis and discarding bad quality data
for XX = 114:150
    biophysics{cellIX}.Rin(biophysics{cellIX}.Rin ==0 ) = NaN;
    biophysics{cellIX}.Rs(biophysics{cellIX}.Rs ==0 ) = NaN;
    biophysics{cellIX}.Cm(biophysics{cellIX}.Cm ==0 ) = NaN;
    biophysics{cellIX}.tau(biophysics{cellIX}.tau ==0 ) = NaN;
end

clear goodXX
for j = 115:150
    [nanmedian(biophysics{j}.Rin(:))/1e9,nanmedian(biophysics{j}.Rs(:)),nanmedian(biophysics{j}.tau(:))*1000,j]
    datasetSingleCells{j}
    goodXX{j} = inputdlg({'Good','TimeConstant'},'UI / Cell selection',1);
end

for j = 115:150
    
    if strcmp(goodXX{j}(1),'0')
        biophysics{j}.Rin = NaN;
    end
    if strcmp(goodXX{j}(2),'0')
        biophysics{j}.tau = NaN;
        biophysics{j}.Cm = NaN;
    end
end
  





%% plot results (pooled from different approaches) that include test pulses




clear Rin_pooled Rs_pooled Cm_pooled tau_pooled
criterion = Anatomicallocationsoverview(:,2) == 1; feature = 1;
% | Anatomicallocationsoverview(:,3) == 1; % only Dp
criterion = Anatomicallocationsoverview(:,3) == 1; feature = 2;
criterion = Anatomicallocationsoverview(:,4) == 1; feature = 3;
criterion = Anatomicallocationsoverview(:,7) == 1; feature = 4;


counter = 1;
for cellIX = 1:150
    if criterion(cellIX)== 1 && ~isempty(biophysics{cellIX})
        Rin_pooled{feature}(counter) = nanmedian(biophysics{cellIX}.Rin(:));
        Rs_pooled{feature}(counter) = nanmedian(biophysics{cellIX}.Rs(:));
        Cm_pooled{feature}(counter) = nanmedian(biophysics{cellIX}.Cm(:));
        tau_pooled{feature}(counter) = nanmedian(biophysics{cellIX}.tau(:));
        counter = counter + 1;
    end
end
counter-1

figure(41); 
subplot(1,3,1); bar([nanmedian(Rin_pooled{1}),nanmedian(Rin_pooled{2}),nanmedian(Rin_pooled{3}),nanmedian(Rin_pooled{4})])
subplot(1,3,2); bar([nanmedian(Rs_pooled{1}),nanmedian(Rs_pooled{2}),nanmedian(Rs_pooled{3}),nanmedian(Rs_pooled{4})])
subplot(1,3,3); bar(1e3*[nanmedian(tau_pooled{1}),nanmedian(tau_pooled{2}),nanmedian(tau_pooled{3}),nanmedian(tau_pooled{4})])

figure(43); 
subplot(1,3,1); hold on; errorbar([nanmedian(Rin_pooled{1}),nanmedian(Rin_pooled{2}),nanmedian(Rin_pooled{3}),nanmedian(Rin_pooled{4})]/1e9,[nanstd(Rin_pooled{1}),nanstd(Rin_pooled{2}),nanstd(Rin_pooled{3}),nanstd(Rin_pooled{4})]/1e9,'.k','Markersize',18,'LineWidth',2)
set(gca,'XTick',[1 2 3 4]); set(gca,'XTickLabel',{'pDp','border','NT','aDp'}); xlim([0.65 4.35])
[sum(~isnan(Rin_pooled{1})),sum(~isnan(Rin_pooled{2})),sum(~isnan(Rin_pooled{3})),sum(~isnan(Rin_pooled{4}))]
subplot(1,3,2); hold on; errorbar([nanmedian(Rs_pooled{1}),nanmedian(Rs_pooled{2}),nanmedian(Rs_pooled{3}),nanmedian(Rs_pooled{4})],[nanmedian(Rs_pooled{1}),nanstd(Rs_pooled{2}),nanstd(Rs_pooled{3}),nanstd(Rs_pooled{4})],'.k','Markersize',18,'LineWidth',2)
set(gca,'XTick',[1 2 3 4]); set(gca,'XTickLabel',{'pDp','border','NT','aDp'}); xlim([0.65 4.35])
[sum(~isnan(Rs_pooled{1})),sum(~isnan(Rs_pooled{2})),sum(~isnan(Rs_pooled{3})),sum(~isnan(Rs_pooled{4}))]
subplot(1,3,3); hold on; errorbar(1e3*[nanmedian(tau_pooled{1}),nanmedian(tau_pooled{2}),nanmedian(tau_pooled{3}),nanmedian(tau_pooled{4})],1e3*[nanstd(tau_pooled{1}),nanstd(tau_pooled{2}),nanstd(tau_pooled{3}),nanstd(tau_pooled{4})],'.k','Markersize',18,'LineWidth',2)
set(gca,'XTick',[1 2 3 4]); set(gca,'XTickLabel',{'pDp','border','NT','aDp'}); xlim([0.65 4.35])
[sum(~isnan(tau_pooled{1})),sum(~isnan(tau_pooled{2})),sum(~isnan(tau_pooled{3})),sum(~isnan(tau_pooled{4}))]




xValues = [ones(size((Rin_pooled{1}))),2*ones(size((Rin_pooled{2}))),3*ones(size((Rin_pooled{3}))),4*ones(size((Rin_pooled{4})))]-0.2+0.4*[rand(size((Rin_pooled{1}))),rand(size((Rin_pooled{2}))),rand(size((Rin_pooled{3}))),rand(size((Rin_pooled{4})))];

figure(43);
subplot(1,3,1); plot(xValues,[(Rin_pooled{1}),(Rin_pooled{2}),(Rin_pooled{3}),(Rin_pooled{4})]/1e9,'.','Color',[0.6 0.6 0.6],'Markersize',9,'LineWidth',2)
set(gca,'XTick',[1 2 3 4]); set(gca,'XTickLabel',{'pDp','border','NT','aDp'}); xlim([0.65 4.35])
subplot(1,3,2); plot(xValues,[(Rs_pooled{1}),(Rs_pooled{2}),(Rs_pooled{3}),(Rs_pooled{4})],'.','Color',[0.6 0.6 0.6],'Markersize',9,'LineWidth',2)
set(gca,'XTick',[1 2 3 4]); set(gca,'XTickLabel',{'pDp','border','NT','aDp'}); xlim([0.65 4.35])
subplot(1,3,3); plot(xValues,1e3*[(tau_pooled{1}),(tau_pooled{2}),(tau_pooled{3}),(tau_pooled{4})],'.','Color',[0.6 0.6 0.6],'Markersize',9,'LineWidth',2)
set(gca,'XTick',[1 2 3 4]); set(gca,'XTickLabel',{'pDp','border','NT','aDp'}); xlim([0.65 4.35])



sample1 = tau_pooled{1}; sample1 = sample1(~isnan(sample1));
sample2 = tau_pooled{4}; sample2 = sample2(~isnan(sample2));


[~,p,~] = kstest2(sample1,sample2);




