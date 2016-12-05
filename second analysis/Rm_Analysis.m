
% AAAB : [5 6 7],[8 9 10 14 15 16],[17 18 19]
% 
% AAAC :  [5 6 7 8],[9 10 11 14 15 16 17]
% 
% AAAD : [4 5 6 7],[8 9 (9)]
% 
% 
% j =3
% % inhTrials = [5 6 7];
% excTrials = [8 9 10 14 15 16];
% excTrials = [9 10 11 14 15 16 17];
% excTrials = [8 9];
% excTrials = [17 18 19]
% 
% FileList = dir('**.xsg');
% clear INH EXC
% counter = 1;
% for i = inhTrials
%     load(FileList(i).name,'-mat');
%     A = data.ephys.trace_1;
%     INH(:,counter) = A;
%     counter = counter + 1;
% end
% counter = 1;
% for i = excTrials
%     load(FileList(i).name,'-mat');
%     A = data.ephys.trace_1;
%     EXC(:,counter) = A;
%     counter = counter + 1;
% end
% 
% 
% excMean = mean(EXC');
%  
% excMean = reshape(excMean,1000,380);
% 
% pixelz = 1;
% baseline = (conv2(excMean(600:990,2:end),fspecial('gaussian',[1 pixelz],3),'same') + conv2(excMean(600:990,1:end-1),fspecial('gaussian',[1 pixelz],3),'same'))/2;
% inline = conv2(excMean(100:490,2:end),fspecial('gaussian',[1 pixelz],3),'same');
% 
% 
% Rin = 5./( median(baseline) - median(inline) );
% 
% 
% IX = 150
% 
% Rs = 45e6;
% 
% biophysics{IX}.Rin = Rin*1e9;
% biophysics{IX}.Rs = Rs;
% 
% Y = mean(excMean');
% X = (0:numel(Y)-1)/10;
% [xData, yData] = prepareCurveData( X(2:500), Y(2:500) );
% ft = fittype( 'a*exp(-x/b)+d*exp(-x/e)+c', 'independent', 'x', 'dependent', 'y' );
% opts = fitoptions( ft ); opts.Display = 'Off';
% opts.Lower = [-Inf -Inf -Inf -Inf -Inf];
% opts.StartPoint = [0.107324056896322 20 0.374643668266154 -22 1];
% opts.Upper = [Inf Inf Inf Inf Inf];
% [fitresult, gof] = fit( xData, yData, ft, opts );
% 
% tauVC = max(fitresult.b,fitresult.e)
% 
% tauVC/Rs*1e9
% 
% biophysics{IX}.Cm = tauVC*1e-3/Rs;
% biophysics{IX}.tau = tauVC*1e-3/Rs*mean(Rin(:))*1e9;


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

            Y = mean(excerpt');
            X = (0:numel(Y)-1)/10;
            [xData, yData] = prepareCurveData( X(2:500), Y(2:500) );
            ft = fittype( 'a*exp(-x/b)+d*exp(-x/e)+c', 'independent', 'x', 'dependent', 'y' );
            opts = fitoptions( ft ); opts.Display = 'Off';
            opts.Lower = [-Inf -Inf -Inf -Inf -Inf];
            opts.StartPoint = [-9 5 -22 -80 0.4];
            opts.Upper = [Inf Inf Inf Inf Inf];
            try
                [fitresult, gof] = fit( xData, yData, ft, opts );

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


for XX = 114:150
    biophysics{cellIX}.Rin(biophysics{cellIX}.Rin ==0 ) = NaN;
    biophysics{cellIX}.Rs(biophysics{cellIX}.Rs ==0 ) = NaN;
    biophysics{cellIX}.Cm(biophysics{cellIX}.Cm ==0 ) = NaN;
    biophysics{cellIX}.tau(biophysics{cellIX}.tau ==0 ) = NaN;
end



clear Rin_pooled Rs_pooled Cm_pooled tau_pooled
criterion = Anatomicallocationsoverview(:,2) == 1;% | Anatomicallocationsoverview(:,3) == 1; % only Dp
criterion = Anatomicallocationsoverview(:,3) == 1;
criterion = Anatomicallocationsoverview(:,4) == 1;
criterion = Anatomicallocationsoverview(:,7) == 1;

feature = 4;

counter = 1;
for cellIX = 115:150
    if criterion(cellIX) == 1
%         Rin_pooled{feature}(counter) = nanmedian(biophysics{cellIX}.Rin(:));
%         Rs_pooled{feature}(counter) = nanmedian(biophysics{cellIX}.Rs(:));
%         Cm_pooled{feature}(counter) = nanmedian(biophysics{cellIX}.Cm(:));
%         tau_pooled{feature}(counter) = nanmedian(biophysics{cellIX}.tau(:));
        counter = counter + 1;
    end
end
counter

figure(41); 
subplot(1,3,1); bar([nanmedian(Rin_pooled{1}),nanmedian(Rin_pooled{2}),nanmedian(Rin_pooled{3}),nanmedian(Rin_pooled{4})])
subplot(1,3,2); bar([nanmedian(Rs_pooled{1}),nanmedian(Rs_pooled{2}),nanmedian(Rs_pooled{3}),nanmedian(Rs_pooled{4})])
subplot(1,3,3); bar(1e3*[nanmedian(tau_pooled{1}),nanmedian(tau_pooled{2}),nanmedian(tau_pooled{3}),nanmedian(tau_pooled{4})])

figure(42); 
subplot(1,3,1); errorbar([nanmedian(Rin_pooled{1}),nanmedian(Rin_pooled{2}),nanmedian(Rin_pooled{3}),nanmedian(Rin_pooled{4})],[nanstd(Rin_pooled{1}),nanstd(Rin_pooled{2}),nanstd(Rin_pooled{3}),nanstd(Rin_pooled{4})],'.k','Markersize',18,'LineWidth',2)
set(gca,'XTick',[1 2 3 4]); set(gca,'XTickLabel',{'pDp','border','NT','aDp'}); xlim([0.75 4.25])
subplot(1,3,2); errorbar([nanmedian(Rs_pooled{1}),nanmedian(Rs_pooled{2}),nanmedian(Rs_pooled{3}),nanmedian(Rs_pooled{4})],[nanmedian(Rs_pooled{1}),nanstd(Rs_pooled{2}),nanstd(Rs_pooled{3}),nanstd(Rs_pooled{4})],'.k','Markersize',18,'LineWidth',2)
set(gca,'XTick',[1 2 3 4]); set(gca,'XTickLabel',{'pDp','border','NT','aDp'}); xlim([0.75 4.25])
subplot(1,3,3); errorbar(1e3*[nanmedian(tau_pooled{1}),nanmedian(tau_pooled{2}),nanmedian(tau_pooled{3}),nanmedian(tau_pooled{4})],1e3*[nanstd(tau_pooled{1}),nanstd(tau_pooled{2}),nanstd(tau_pooled{3}),nanstd(tau_pooled{4})],'.k','Markersize',18,'LineWidth',2)
set(gca,'XTick',[1 2 3 4]); set(gca,'XTickLabel',{'pDp','border','NT','aDp'}); xlim([0.75 4.25])









