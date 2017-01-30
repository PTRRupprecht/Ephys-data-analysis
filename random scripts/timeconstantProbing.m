

for j = 1:7

    if j == 1
        cd('C:\Data\rupppete\PhD\electrophysiology2016\SingleCells\161201_AAAB');
        excTrials = [8 9 10 14 15 16];
    elseif j == 2
        cd('C:\Data\rupppete\PhD\electrophysiology2016\SingleCells\161201_AAAC');
        excTrials = [9 10 11 14 15 16 17];
    elseif j == 3
        cd('C:\Data\rupppete\PhD\electrophysiology2016\SingleCells\161201_AAAD');
        excTrials = [8 9 ];
    elseif j == 4
        cd('C:\Data\rupppete\PhD\electrophysiology2016\SingleCells\161211_AAAA');
%         excTrials = [ 7 8 9 (10) (11) 12 13 (14) 15];% 15 16 17 18 ];
        excTrials = [ 7 8 9 12 13 15];% 15 16 17 18 ];
    elseif j == 5
        cd('C:\Data\rupppete\PhD\electrophysiology2016\SingleCells\161211_AAAB');
%         excTrials = [7 8 (9)];% 15 16 17 18 ];
        excTrials = [7 8];% 15 16 17 18 ];
    elseif j == 6
        cd('C:\Data\rupppete\PhD\electrophysiology2016\SingleCells\161211_AAAC');
%         excTrials = [4 5 ((6)) ((7))];% 15 16 17 18 ];
        excTrials = [4 5 ];% 15 16 17 18 ];
    elseif j == 7
        cd('C:\Data\rupppete\PhD\electrophysiology2016\SingleCells\161211_AAAD');
%         excTrials = [4 5 ((6))];% 15 16 17 18 ];
        excTrials = [4 5];% 15 16 17 18 ];
    end

    FileList = dir('**.xsg');
    clear EXC
    counter = 1;
    for i = excTrials
        load(FileList(i).name,'-mat');
        A = data.ephys.trace_1;
        EXC(:,counter) = A;
        counter = counter + 1;
    end
    excMean0 = mean(EXC');
    for k = 1:size(excMean0,2)/1000
        k
        
        excMean = excMean0;
        excMean = reshape(excMean,1000,380);
        try
            excMeanpp = excMean(:,k+1);
        catch
            excMeanpp = excMean(:,k);
        end
        try
            excMeanppp = excMean(:,k+2);
        catch
            excMeanppp = excMeanpp;
        end
        try
            excMeanppp = excMean(:,k+2);
        catch
            excMeanppp = excMeanpp;
        end
        excMean = excMean(:,k);
        pixelz = 1;
        baseline = (conv2(excMean(600:990,1),fspecial('gaussian',[1 pixelz],3),'same') + conv2(excMeanpp(600:990,1),fspecial('gaussian',[1 pixelz],3),'same') + conv2(excMeanppp(600:990,1),fspecial('gaussian',[1 pixelz],3),'same'))/3;
        inline = (conv2(excMean(100:490,1),fspecial('gaussian',[1 pixelz],3),'same') + conv2(excMeanpp(100:490,1),fspecial('gaussian',[1 pixelz],3),'same'))/2;

        Rin = 5./( median(baseline) - median(inline) );
        RinX(j,k) = Rin;
        Rs1(j,k) = 5e-3./(median(excMean)-min(excMean))*1e12/1e6;
        Rs2(j,k) = - 5e-3./(median(excMean)-max(excMean))*1e12/1e6;

        %     biophysics{IX}.Rin = Rin*1e9;
        %     biophysics{IX}.Rs = [Rs1,Rs2];

        Y = excMean;
        X = (0:numel(Y)-1)/10;
%         diffX = Y - circshift(Y,[-1]); startIDX = find(diffX<0,1,'first');
%         [xData, yData] = prepareCurveData( X(startIDX:500), Y(startIDX:500)' );
%         ft = fittype( 'a*exp(-x/b)+d*exp(-x/e)+c', 'independent', 'x', 'dependent', 'y' );
%         opts = fitoptions( ft ); opts.Display = 'Off';
%         opts.Lower = [-Inf -Inf -Inf -Inf -Inf];
%         opts.StartPoint = [-20 4 -20 -22 0.1];
%         opts.Upper = [Inf Inf Inf Inf Inf];
%         try
%             [fitresult, gof] = fit( xData, yData, ft, opts );
%             tauX(j,k) = max(fitresult.b,fitresult.e)/median(median([Rs1(j,k) Rs2(j,k)]).*Rin*1e9/1e6);
%         catch
%             tauX(j,k) = NaN;
%         end
    end

    %     biophysics{IX}.Cm = tauVC*1e-3/mean([Rs1 Rs2]);
    %     biophysics{IX}.tau = tauVC*1e-3/mean([Rs1 Rs2])*mean(Rin(:))*1e9;
    
end

for j = 1:7

    if j == 1
        cd('C:\Data\rupppete\PhD\electrophysiology2016\SingleCells\161201_AAAB');
        excTrials = [5 6 7];
    elseif j == 2
        cd('C:\Data\rupppete\PhD\electrophysiology2016\SingleCells\161201_AAAC');
        excTrials = [5 6 7 8];
    elseif j == 3
        cd('C:\Data\rupppete\PhD\electrophysiology2016\SingleCells\161201_AAAD');
        excTrials = [4 5 6 7];
    elseif j == 4
        cd('C:\Data\rupppete\PhD\electrophysiology2016\SingleCells\161211_AAAA');
%         excTrials = [(3) 4 5];% 15 16 17 18 ];
        excTrials = [4 5];% 15 16 17 18 ];
    elseif j == 5
        cd('C:\Data\rupppete\PhD\electrophysiology2016\SingleCells\161211_AAAB');
%         excTrials = [3 4 (5) 6];% 15 16 17 18 ];
        excTrials = [3 4 6];% 15 16 17 18 ];
    elseif j == 6
        cd('C:\Data\rupppete\PhD\electrophysiology2016\SingleCells\161211_AAAC');
        excTrials = [2 3];% 15 16 17 18 ];
    elseif j == 7
        cd('C:\Data\rupppete\PhD\electrophysiology2016\SingleCells\161211_AAAD');
        excTrials = [2 3];% 15 16 17 18 ];
    end

    FileList = dir('**.xsg');
    clear INH
    counter = 1;
    for i = excTrials
        load(FileList(i).name,'-mat');
        A = data.ephys.trace_1;
        INH(:,counter) = A;
        counter = counter + 1;
    end
    inhMean0(j,:) = mean(INH');
end



cmap = distinguishable_colors(7);
for j = 1:7

    if j == 1
        cd('C:\Data\rupppete\PhD\electrophysiology2016\SingleCells\161201_AAAB');
        excTrials = [8 9 10 14 15 16];
    elseif j == 2
        cd('C:\Data\rupppete\PhD\electrophysiology2016\SingleCells\161201_AAAC');
        excTrials = [9 10 11 14 15 16 17];
    elseif j == 3
        cd('C:\Data\rupppete\PhD\electrophysiology2016\SingleCells\161201_AAAD');
        excTrials = [8 9 ];
    elseif j == 4
        cd('C:\Data\rupppete\PhD\electrophysiology2016\SingleCells\161211_AAAA');
%         excTrials = [ 7 8 9 (10) (11) 12 13 (14) 15];% 15 16 17 18 ];
        excTrials = [ 7 8 9 12 13 15];% 15 16 17 18 ];
    elseif j == 5
        cd('C:\Data\rupppete\PhD\electrophysiology2016\SingleCells\161211_AAAB');
%         excTrials = [7 8 (9)];% 15 16 17 18 ];
        excTrials = [7 8];% 15 16 17 18 ];
    elseif j == 6
        cd('C:\Data\rupppete\PhD\electrophysiology2016\SingleCells\161211_AAAC');
%         excTrials = [4 5 ((6)) ((7))];% 15 16 17 18 ];
        excTrials = [4 5 ];% 15 16 17 18 ];
    elseif j == 7
        cd('C:\Data\rupppete\PhD\electrophysiology2016\SingleCells\161211_AAAD');
%         excTrials = [4 5 ((6))];% 15 16 17 18 ];
        excTrials = [4 5];% 15 16 17 18 ];
    end

    FileList = dir('**.xsg');
    clear EXC
    counter = 1;
    for i = excTrials
        load(FileList(i).name,'-mat');
        A = data.ephys.trace_1;
        EXC(:,counter) = A;
        counter = counter + 1;
    end
    excMean0 = mean(EXC');
%     figure(6); subplot(3,1,1); hold on; plot((1:numel(excMean0))/1e4,excMean0+j*200,'Color',cmap(j,:))
    
    excMean = reshape(excMean0,1000,380);
    excMean = nanmedian(excMean,2);
    excMean = repmat(excMean,[1 380]);
    excMean = reshape(excMean,[1 380000]);
    figure(6); subplot(3,1,2); hold on; plot((1:numel(excMean0))/1e4,excMean0-excMean+j*43,'Color',cmap(j,:))
    figure(6); subplot(3,1,2); hold on; plot((1:numel(inhMean0(j,:)))/1e4,inhMean0(j,:)+j*43+350,'Color',cmap(j,:))
    figure(6); subplot(3,1,1); plot((1:numel(inhMean0(j,:)))/1e4,1./smooth(1./median(RinX(j,:))+((inhMean0(j,:)-median(inhMean0(j,:)))-(excMean0-excMean))/70+j,100),'Color',cmap(j,:)); hold on;
    ylim([0 2])
end

% for j == 1 !
t1 = mean(excMean(:,119:124),2);
t1 = mean(excMean(:,121:125),2);
t2 = mean(excMean(:,1:100),2);
1/([mean(t1(730:999))-mean(t1(230:499))]*1e-12/5e-3)
1/([mean(t2(730:999))-mean(t2(230:499))]*1e-12/5e-3)
tix = (1:1000)/10;
tii = (-500:1:-1)/10;
figure; plot([tii,tix],[smooth(t1(501:end))- mean(t1(730:999));smooth(t1,10) - mean(t1(730:999))],'r'),
hold on;
plot([tii,tix],[smooth(t1(501:end))- mean(t1(730:999));smooth(t2,10) - mean(t2(730:999))],'k'),


cmap = distinguishable_colors(7);
for j = [1 2 3 4 5 6 7]
%     if j ==5
%         tauXX = nanmedian(tauX(1:3,:),1);
%     else
%         tauXX = tauX(j,:);
%     end
%     tauZ = ([tauXX; circshift(tauXX,[0 1]); circshift(tauXX,[0 2]); circshift(tauXX,[0 3]); circshift(tauXX,[0 4]); circshift(tauXX,[0 5]); circshift(tauXX,[0 6]); circshift(tauXX,[0 7]); circshift(tauXX,[0 8])]);
%     timet = (1:size(tauZ,2))/10-0.5;
%     figure(6), subplot(3,1,3); hold on; plot(timet,nanmedian(tauZ),'.','Color',cmap(j,:)); ylim([0 500])
    if 0%j ==5
        tauXX = nanmedian(tauX(1:3,:),1);
    else
        RinXX = RinX(j,:);
    end
    RinZ = ([RinXX; circshift(RinXX,[0 1]); circshift(RinXX,[0 2]); circshift(RinXX,[0 3]); circshift(RinXX,[0 4])]);% circshift(RinXX,[0 5]); circshift(RinXX,[0 6]); circshift(RinXX,[0 7]); circshift(RinXX,[0 8])]);
    if j < 4
        timet = (1:size(RinZ,2))/10 -0.2;
    else
        timet = (1:size(RinZ,2))/10 -0.2 +0.2;
    end
    RinZnn = nanmedian(RinZ(:));
    figure(991), subplot(3,1,3); hold on; plot(timet,nanmedian(RinZ)/RinZnn,'.','Color',cmap(j,:)); ylim([0 2])
end


clear RinDuringDis RinBeforeDis
for j = 1:7
    RinBefore(j) = median(RinX(j,:));
    RinBeforeStd(j) = std(RinX(j,:));
    RinBeforeDis(j,:) = (RinX(j,:));
    if j < 4
        RinDuring(j) = median(RinX(j,120:124));
        RinDuringStd(j) = std(RinX(j,120:124));
        RinDuringDis(j,:) = (RinX(j,120:124));
    else
        RinDuring(j) = median(RinX(j,117:121));
        RinDuringStd(j) = std(RinX(j,117:121));
        RinDuringDis(j,:) = (RinX(j,117:121));
    end
end


figure(5421), bar([1./RinBefore([1:3 5:7])',1./RinDuring([1:3 5:7])'])
figure(541), bar([RinBefore([1:3 5:7])',RinDuring([1:3 5:7])'])
figure(5431), bar([mean(RinBefore([1:3 5:7])), mean(RinDuring([1:3 5:7])')])
xlim([0.25 2.75]);
hold on;
for j = [1 2 3 5 6 7]
    plot([1 2],[mean(RinBefore(j)), mean(RinDuring(j)')],'-')
end