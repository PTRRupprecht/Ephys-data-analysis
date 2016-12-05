%% quantify area-under-curve after supposed response onset in a 1.5 sec-time window
% this script uses variables created by "TimingAnalysis.m"
% Basically, the script measures the difference in the AUC between trials
% of the same neuron with the same odor vs. trials with different odors.
% The measurement of the difference is done via absolute difference
% (subtraction).
% Afterwards, comparisons of those differences within odors (=variability,
% 'intra') with those across odors can be made, all within one neuron.

counter = 1;
for IX = 1:98
    IX
    
    if 1 % answer{IX}{1}
        onsetVC70 = X_all{IX}.onsetVC70; onsetVC70(onsetVC70 == 0) = NaN;
        onsetVC00 = X_all{IX}.onsetVC00; onsetVC00(onsetVC00 == 0) = NaN;
        odorsVC = X_all{IX}.odorsVC;

        for p = 1:numel(odorsVC)
            odorIndex = find(~cellfun(@isempty,strfind(datasetSingleCells{IX}.odors,odorsVC{p})));
            dayIndex = find(~cellfun(@isempty,strfind(days,X_all{IX}.dateID)));

            trials_T = find(~cellfun(@isempty,strfind(datasetSingleCells{IX}.VC0odor,X_all{IX}.odorsVC{p}) ));
            trials = datasetSingleCells{IX}.VC0(trials_T);
            GoodTrials = squeeze(chosen_ones(IX,:,:));
            takeLeave00 = [];
            for jj = 1:numel(trials)
                takeLeave00(jj) = any(GoodTrials(:)==trials(jj));
            end
            trials_T = find(~cellfun(@isempty,strfind(datasetSingleCells{IX}.VC70odor,X_all{IX}.odorsVC{p}) ));
            trials = datasetSingleCells{IX}.VC70(trials_T);
            GoodTrials = squeeze(chosen_ones(IX,:,:));
            takeLeave70 = [];
            for jj = 1:numel(trials)
                takeLeave70(jj) = any(GoodTrials(:)==trials(jj));
            end
            takeLeave70 = find(~takeLeave70);
            takeLeave00 = find(~takeLeave00);
            onsetVC70(p,takeLeave70,:) = NaN;
            onsetVC00(p,takeLeave00,:) = NaN;

            delay = onsetZ(dayIndex,odorIndex);
            onsetVC70(p,:,:) = circshift(onsetVC70(p,:,:),[0 0 1e5-round(delay)]);
            try
                onsetVC00(p,:,:) = circshift(onsetVC00(p,:,:),[0 0 1e5-round(delay)]);
            end
        end




        intraDiff_inh{IX} = []; 
        intraDiff_exc{IX} = [];
        interDiff_exc{IX} = [];
        interDiff_inh{IX} = []; 
        strengthIntraI{IX} = [];
        strengthInterI{IX} = [];
        strengthIntraE{IX} = [];
        strengthInterE{IX} = [];

        firstIndex = repmat([1:size(onsetVC70,1)]',[1 size(onsetVC70,2)]);
        nanIndex = ~isnan(squeeze(mean(onsetVC70,3)));
        for kk = 1:(numel(firstIndex)-1)

            for dd = (kk+1):numel(firstIndex)
                if nanIndex(kk) && nanIndex(dd)
                    p1 = firstIndex(kk);
                    p2 = firstIndex(dd);

                    for tt = 40
                        t_end = 1e5+tt*500-500;
                        t_start = 1e5-500; % IMPORTANT : either integrate or moving window

                        [x,y] = ind2sub([size(firstIndex,1) size(firstIndex,2)],kk);
                        [x2,y2] = ind2sub([size(firstIndex,1) size(firstIndex,2)],dd);

                        y2 = find(nanIndex(x2,:));
                        if p1 == p2; y2 = setdiff(y2,y); end

                        AUCE1 = nanmean(nanmean(squeeze(onsetVC70(x,y,t_start:t_end))) - squeeze(nanmean(onsetVC70(x,y,1:0.9e5),3)) );
                        AUCE2 = nanmean( nanmean(squeeze(onsetVC70(x2,y2,t_start:t_end)),2) - squeeze(nanmean(onsetVC70(x2,y2,1:0.9e5),3))' );

                        if p1 == p2 % intra
                            intraDiff_exc{IX} = [intraDiff_exc{IX}, abs(AUCE1-AUCE2)];
                            strengthIntraE{IX} = [strengthIntraE{IX}, abs(AUCE1+AUCE2)];
                        else % inter
                            interDiff_exc{IX} = [interDiff_exc{IX}, abs(AUCE1-AUCE2)];
                            strengthInterE{IX} = [strengthInterE{IX}, abs(AUCE1+AUCE2)];
                        end


                    end

                end
            end
        end

        firstIndex = repmat([1:size(onsetVC00,1)]',[1 size(onsetVC00,2)]);
        nanIndex = ~isnan(squeeze(mean(onsetVC00,3)));
        for kk = 1:(numel(firstIndex)-1)

            for dd = (kk+1):numel(firstIndex)
                if nanIndex(kk) && nanIndex(dd)
                    p1 = firstIndex(kk);
                    p2 = firstIndex(dd);

                    for tt = 30
                        t_end = 1e5+tt*500;
                        t_start = 1e5-500; % IMPORTANT : either integrate or moving window

                        [x,y] = ind2sub([size(firstIndex,1) size(firstIndex,2)],kk);
                        [x2,y2] = ind2sub([size(firstIndex,1) size(firstIndex,2)],dd);

                        y2 = find(nanIndex(x2,:));
                        if p1 == p2; y2 = setdiff(y2,y); end

                        try
                            AUCI1 = nanmean(nanmean(squeeze(onsetVC00(x,y,t_start:t_end))) - squeeze(nanmean(onsetVC00(x,y,1:0.9e5),3)) );
                            AUCI2 = nanmean(nanmean(squeeze(onsetVC00(x2,y2,t_start:t_end)),2) - squeeze(nanmean(onsetVC00(x2,y2,1:0.9e5),3)' ) );
                        catch
                            disp(['Problem with cell no.',32,num2str(IX)]);
                        end

                        if p1 == p2 % intra
                            intraDiff_inh{IX} = [intraDiff_inh{IX}, abs(AUCI1-AUCI2)];
                            strengthIntraI{IX} = [strengthIntraI{IX}, abs(AUCI1+AUCI2)];
                        else % inter
                            interDiff_inh{IX} = [interDiff_inh{IX}, abs(AUCI1-AUCI2)];
                            strengthInterI{IX} = [strengthInterI{IX}, abs(AUCI1+AUCI2)];
                        end


                    end

                end
            end
        end

    end
    
    
end


%% plot data points and distribution
figure(412); 

XX = []; YY = []; XX2 = []; YY2 = [];
rng(77);
for IX = 1:98
    if  answer_2{IX} == 2 %str2num(answer{IX}{1}) == 1
        XX(IX) = mean(interDiff_inh{IX}) - mean(intraDiff_inh{IX});
        WW(IX) = sqrt(numel(strengthIntraI{IX})+numel(strengthInterI{IX}));%*(mean(strengthIntraI{IX})+mean(strengthInterI{IX}))/2/9;
        XX2(IX) = mean(interDiff_exc{IX}) - mean(intraDiff_exc{IX});
        WW2(IX) = sqrt(numel(strengthIntraE{IX})+numel(strengthInterE{IX}));%*(mean(strengthIntraE{IX})+mean(strengthInterE{IX}))/2/9;
    %     plot(XX(IX),'.','MarkerSize',max(2,WW(IX))); hold on;
        offset1 = rand()*0.03;
    %     if WW(IX) > 10 || WW2(IX) > 10
    %     plot([XX(IX),XX2(IX)],[0.25+offset1,0.3+offset1],'Color',[0.8 0.8 0.8],'LineWidth',sqrt(max(0,WW2(IX))*max(0,WW(IX)))/10 ); hold on;
    %     end
        plot(XX(IX),0.25+offset1,'.r','MarkerSize',3*sqrt(max(2,WW(IX))*max(2,WW(IX)))); hold on;
        plot(XX2(IX),0.30+offset1,'.k','MarkerSize',3*sqrt(max(2,WW2(IX))*max(2,WW2(IX)))); hold on;

    %     plot(0:0.1:20,0:0.1:20,'k');
    end
end
WW(isnan(XX)) = [];
XX(isnan(XX)) = [];
FXX = ksdensity(XX,-8:0.01:12,'weights',WW.^2,'bandwidth',0.9);
WW2(isnan(XX2)) = [];
XX2(isnan(XX2)) = [];
FXX2 = ksdensity(XX2,-8:0.01:12,'weights',WW2.^2,'bandwidth',0.9);
plot(0*(-9:1:10),(0:19)/50,'k')
plot(-8:0.01:12,FXX,'r'); hold on; plot(-8:0.01:12,FXX2,'k'); hold off;


find(XX> -2.249 & XX < -2.247)
manualSelection(38,X_all,datasetSingleCells,days,onsetZ)




%% variability explained paradigm
counter = 1;
clear AUC_inh AUC_exc AUC_odor
clear trialsInh trialsExc
for IX = 1:98
    
        onsetVC70 = X_all{IX}.onsetVC70; onsetVC70(onsetVC70 == 0) = NaN;
        onsetVC00 = X_all{IX}.onsetVC00; onsetVC00(onsetVC00 == 0) = NaN;
        odorsVC = X_all{IX}.odorsVC;

        for p = 1:numel(odorsVC)
            odorIndex = find(~cellfun(@isempty,strfind(datasetSingleCells{IX}.odors,odorsVC{p})));
            dayIndex = find(~cellfun(@isempty,strfind(days,X_all{IX}.dateID)));

            trials_T = find(~cellfun(@isempty,strfind(datasetSingleCells{IX}.VC0odor,X_all{IX}.odorsVC{p}) ));
            trials = datasetSingleCells{IX}.VC0(trials_T);
            if ~isempty(trials)
                trialsInh{IX}(p,1:numel(trials)) = trials;
            end
            GoodTrials = squeeze(chosen_ones(IX,:,:));
            takeLeave00 = [];
            for jj = 1:numel(trials)
                takeLeave00(jj) = any(GoodTrials(:)==trials(jj));
            end
            trials_T = find(~cellfun(@isempty,strfind(datasetSingleCells{IX}.VC70odor,X_all{IX}.odorsVC{p}) ));
            trials = datasetSingleCells{IX}.VC70(trials_T);
            trialsExc{IX}(p,1:numel(trials)) = trials;
            GoodTrials = squeeze(chosen_ones(IX,:,:));
            takeLeave70 = [];
            for jj = 1:numel(trials)
                takeLeave70(jj) = any(GoodTrials(:)==trials(jj));
            end
            takeLeave70 = find(~takeLeave70);
            takeLeave00 = find(~takeLeave00);
            onsetVC70(p,takeLeave70,:) = NaN;
            onsetVC00(p,takeLeave00,:) = NaN;

            delay = onsetZ(dayIndex,odorIndex);
            onsetVC70(p,:,:) = circshift(onsetVC70(p,:,:),[0 0 1e5-round(delay)]);
            try
                onsetVC00(p,:,:) = circshift(onsetVC00(p,:,:),[0 0 1e5-round(delay)]);
            end
        end
        t_end = 1e5+40*500-500;
        t_start = 1e5-500; % IMPORTANT : either integrate or moving window
        if ~isempty(onsetVC00)
            for p = 1:size(onsetVC00,1)
                AUC_inh{IX}(p,:) = mean(onsetVC00(p,:,t_start:t_end),3) - mean(onsetVC00(p,:,1:0.9e5),3); 
            end
        else
             AUC_inh{IX}(1,:) = [NaN];
        end
        for p = 1:size(onsetVC70,1)
            AUC_exc{IX}(p,:) = mean(onsetVC70(p,:,t_start:t_end),3) - mean(onsetVC70(p,:,1:0.9e5),3);
            AUC_odor{IX}{p} = odorsVC{p};
        end
    
end


figure(45731)
counter = 1;
clear AUC_excL AUC_inhL var_expl_tuning_inh var_expl_tuning_inh2 var_expl_tuning_exc var_expl_tuning_exc2
for IX = 1:98
    for p = 1:size(AUC_exc{IX},1)
        AUC_excX{IX}(p,:) = AUC_exc{IX}(p,:) - nanmean(AUC_exc{IX}(p,:) ); %#ok<SAGROW>
    end
    for p = 1:size(AUC_inh{IX},1)
        AUC_inhX{IX}(p,:) = AUC_inh{IX}(p,:) - nanmean(AUC_inh{IX}(p,:) ); %#ok<SAGROW>
    end
        
    if  answer_2{IX} == 2 && str2num(answer{IX}{1}) == 1 && size(AUC_exc{IX},1)>1 && sum(sum(~isnan(AUC_exc{IX}),2)>1) == size(AUC_exc{IX},1) %#ok<ST2NM>
        ft = fittype( 'a*x+b', 'independent', 'x', 'dependent', 'y' );
        opts = fitoptions( ft );
        opts.Display = 'Off';

        x = trialsExc{IX}(:);
        y = AUC_exc{IX}(:);
        [xData, yData] = prepareCurveData( x, y );
        [fitresult, gof] = fit( xData, yData, ft, opts );
        AUC_excL{IX} = AUC_exc{IX}(:) - fitresult(x);
    end
    if  answer_2{IX} == 2 && str2num(answer{IX}{1}) == 1 && size(AUC_inh{IX},1)>1 && sum(sum(~isnan(AUC_inh{IX}),2)>1) == size(AUC_inh{IX},1) %#ok<ST2NM>

        x = trialsInh{IX}(:);
        y = AUC_inh{IX}(:);
        [xData, yData] = prepareCurveData( x, y );
        [fitresult, gof] = fit( xData, yData, ft, opts );
        AUC_inhL{IX} = AUC_inh{IX}(:) - fitresult(x);

    end
    
    figure(45731);
    if  answer_2{IX} == 2 && str2num(answer{IX}{1}) == 1 && size(AUC_inh{IX},1)>1 && sum(sum(~isnan(AUC_inh{IX}),2)>1) == size(AUC_inh{IX},1) %#ok<ST2NM>
        
        var_expl_tuning_inh(IX) = (nanvar(AUC_inh{IX}(:)) - nanvar(AUC_inhX{IX}(:)));%/nanvar(AUC_inh{IX}(:)); %#ok<SAGROW>
        var_expl_tuning_inh2(IX) = (nanvar(AUC_inh{IX}(:)) - nanvar(AUC_inhL{IX}(:)));%/nanvar(AUC_inh{IX}(:)); %#ok<SAGROW>
        plot(IX,var_expl_tuning_inh(IX)*100,'.' ); hold on;
        plot(IX,var_expl_tuning_inh2(IX)*100,'x' ); hold on;
    end
    if  answer_2{IX} == 2 && str2num(answer{IX}{1}) == 1 && size(AUC_exc{IX},1)>1 && sum(sum(~isnan(AUC_exc{IX}),2)>1) == size(AUC_exc{IX},1) %#ok<ST2NM>
        
        var_expl_tuning_exc(IX) = (nanvar(AUC_exc{IX}(:)) - nanvar(AUC_excX{IX}(:)));%/nanvar(AUC_exc{IX}(:)); %#ok<SAGROW>
        var_expl_tuning_exc2(IX) = (nanvar(AUC_exc{IX}(:)) - nanvar(AUC_excL{IX}(:)));%/nanvar(AUC_exc{IX}(:)); %#ok<SAGROW>
        plot(IX,var_expl_tuning_exc2(IX)*100,'.k'  ); hold on;
    end
end


% resampling and scrambling
clear var_expl_tuning_excRand  var_expl_tuning_exc2Rand var_expl_tuning_inh2Rand var_expl_tuning_inhRand
for repp = 1:10
    repp
    warning('off')
    counter = 1;
    for IX = 1:98
        % permutations, random
        ixf = randperm(size(AUC_exc{IX}(:)));
        ixf = reshape(ixf,size(AUC_exc{IX}));
        AUC_excRand{IX} = AUC_exc{IX}(ixf);
        ixf = randperm(size(AUC_inh{IX}(:)));
        ixf = reshape(ixf,size(AUC_inh{IX}));
        AUC_inhRand{IX} = AUC_inh{IX}(ixf);
        for p = 1:size(AUC_exc{IX},1)
            AUC_excX{IX}(p,:) = AUC_excRand{IX}(p,:) - nanmean(AUC_excRand{IX}(p,:) ); %#ok<SAGROW>
        end
        for p = 1:size(AUC_inhRand{IX},1)
            AUC_inhX{IX}(p,:) = AUC_inhRand{IX}(p,:) - nanmean(AUC_inhRand{IX}(p,:) ); %#ok<SAGROW>
        end
        if  str2num(answer{IX}{1}) == 1 && size(AUC_excRand{IX},1)>1 && sum(sum(~isnan(AUC_excRand{IX}),2)>1) == size(AUC_excRand{IX},1) %#ok<ST2NM>
            ft = fittype( 'a*x+b', 'independent', 'x', 'dependent', 'y' );
            opts = fitoptions( ft );
            opts.Display = 'Off';
            x = trialsExc{IX}(:);
            y = AUC_excRand{IX}(:);
            [xData, yData] = prepareCurveData( x, y );
            [fitresult, gof] = fit( xData, yData, ft, opts );
            AUC_excL{IX} = AUC_excRand{IX}(:) - fitresult(x);
        end
        if  str2num(answer{IX}{1}) == 1 && size(AUC_inhRand{IX},1)>1 && sum(sum(~isnan(AUC_inhRand{IX}),2)>1) == size(AUC_inhRand{IX},1) %#ok<ST2NM>
            x = trialsInh{IX}(:);
            y = AUC_inhRand{IX}(:);
            [xData, yData] = prepareCurveData( x, y );
            [fitresult, gof] = fit( xData, yData, ft, opts );
            AUC_inhL{IX} = AUC_inhRand{IX}(:) - fitresult(x);
        end
        if  str2num(answer{IX}{1}) == 1 && size(AUC_inhRand{IX},1)>1 && sum(sum(~isnan(AUC_inhRand{IX}),2)>1) == size(AUC_inhRand{IX},1) %#ok<ST2NM>
            var_expl_tuning_inhRand(IX,repp) = (nanvar(AUC_inhRand{IX}(:)) - nanvar(AUC_inhX{IX}(:)));nanvar(AUC_inhRand{IX}(:)); %#ok<SAGROW>
            var_expl_tuning_inh2Rand(IX,repp) = (nanvar(AUC_inhRand{IX}(:)) - nanvar(AUC_inhL{IX}(:)));nanvar(AUC_inhRand{IX}(:)); %#ok<SAGROW>
        end
        if  str2num(answer{IX}{1}) == 1 && size(AUC_excRand{IX},1)>1 && sum(sum(~isnan(AUC_excRand{IX}),2)>1) == size(AUC_excRand{IX},1) %#ok<ST2NM>
            var_expl_tuning_excRand(IX,repp) = (nanvar(AUC_excRand{IX}(:)) - nanvar(AUC_excX{IX}(:)));nanvar(AUC_excRand{IX}(:)); %#ok<SAGROW>
            var_expl_tuning_exc2Rand(IX,repp) = (nanvar(AUC_excRand{IX}(:)) - nanvar(AUC_excL{IX}(:)));nanvar(AUC_excRand{IX}(:)); %#ok<SAGROW>
        end
    end
end

LL = [];
LLix = [];
LL = [LL,var_expl_tuning_inh(var_expl_tuning_inh~=0)];
LLix = [LLix,1*ones(size(var_expl_tuning_inh(var_expl_tuning_inh~=0)))];
LL = [LL,var_expl_tuning_inhRand(var_expl_tuning_inhRand~=0)'];
LLix = [LLix,2*ones(size(var_expl_tuning_inhRand(var_expl_tuning_inhRand~=0)))'];
LL = [LL,var_expl_tuning_exc(var_expl_tuning_exc~=0)];
LLix = [LLix,3*ones(size(var_expl_tuning_exc(var_expl_tuning_exc~=0)))];
LL = [LL,var_expl_tuning_excRand(var_expl_tuning_excRand~=0)'];
LLix = [LLix,4*ones(size(var_expl_tuning_excRand(var_expl_tuning_excRand~=0)))'];
LL = [LL,var_expl_tuning_inh2(var_expl_tuning_inh2~=0)];
LLix = [LLix,5*ones(size(var_expl_tuning_inh2(var_expl_tuning_inh2~=0)))];
LL = [LL,var_expl_tuning_inh2Rand(var_expl_tuning_inh2Rand~=0)'];
LLix = [LLix,6*ones(size(var_expl_tuning_inh2Rand(var_expl_tuning_inh2Rand~=0)))'];
LL = [LL,var_expl_tuning_exc2(var_expl_tuning_exc2~=0)];
LLix = [LLix,7*ones(size(var_expl_tuning_exc2(var_expl_tuning_exc2~=0)))];
LL = [LL,var_expl_tuning_exc2Rand(var_expl_tuning_exc2Rand~=0)'];
LLix = [LLix,8*ones(size(var_expl_tuning_exc2Rand(var_expl_tuning_exc2Rand~=0)))'];

figure, boxplot(100*(LL),LLix,'colorgroup',[1 2 1 2 1 2 1 2],'notch','on','outliersize',1,'labels',{'Tuning inh','','Tuning exc','','Adaptation inh','','Adaptation exc',''});
figure, boxplot(sqrt(LL),LLix,'colorgroup',[1 2 1 2 ],'notch','on','outliersize',1,'labels',{'Tuning inh','','Tuning exc',''});



% statistics of used AUCs
EY = []; for k = 1:98; if answer_2{k} == 2; EY = [EY,nanmean(AUC_inh{k}(:))]; end; end;
EX = []; for k = 1:98; if answer_2{k} == 2; EX = [EX,nanmean(AUC_exc{k}(:))]; end; end;
nanmedian(EY),nanstd(EY)
nanmedian(EX),nanstd(EX)
















figure, subplot(2,2,1); boxplot(var_expl_tuning_inh(var_expl_tuning_inh~=0),1)
subplot(2,2,2); boxplot(var_expl_tuning_exc(var_expl_tuning_exc~=0),1)
subplot(2,2,3); boxplot(var_expl_tuning_inh2(var_expl_tuning_inh2~=0),1)
subplot(2,2,4); boxplot(var_expl_tuning_exc2(var_expl_tuning_exc2~=0),1)

mean(var_expl_tuning_inh(var_expl_tuning_inh~=0)),std(var_expl_tuning_inh(var_expl_tuning_inh~=0))
mean(var_expl_tuning_exc(var_expl_tuning_exc~=0)),std(var_expl_tuning_exc(var_expl_tuning_exc~=0))
mean(var_expl_tuning_inh2(var_expl_tuning_inh2~=0)),std(var_expl_tuning_inh2(var_expl_tuning_inh2~=0))
mean(var_expl_tuning_exc2(var_expl_tuning_exc2~=0)),std(var_expl_tuning_exc2(var_expl_tuning_exc2~=0))
hold off;



%% scatter-plot of inhibitory vs. excitatory tuning
figure(411); 

XX = []; YY = []; XX2 = []; YY2 = []; 
for IX = 1:98
    XX(IX) = mean(interDiff_inh{IX}) - mean(intraDiff_inh{IX});
    WW(IX) = sqrt(numel(strengthIntraI{IX})+numel(strengthInterI{IX}));%*(mean(strengthIntraI{IX})+mean(strengthInterI{IX}))/2/9;
    XX2(IX) = mean(interDiff_exc{IX}) - mean(intraDiff_exc{IX});
    WW2(IX) = sqrt(numel(strengthIntraE{IX})+numel(strengthInterE{IX}));%*(mean(strengthIntraE{IX})+mean(strengthInterE{IX}))/2/9;
    subplot(1,2,1); plot(XX(IX),XX2(IX),'.k','MarkerSize',3*sqrt(max(2,WW(IX))*max(2,WW(IX)))); hold on;
    subplot(1,2,2); plot(XX(IX),XX2(IX),'.k','MarkerSize',15); hold on;
end
hold off;


%% alternative distribution plot based on Anne Urai's violin plot

colors = cbrewer('qual', 'Set1', 10);

XX = []; YY = []; XX2 = []; YY2 = []; 
L1 = [];
L2 = [];
for IX = 1:98
    XX(IX) = mean(interDiff_inh{IX}) - mean(intraDiff_inh{IX});
    WW(IX) = sqrt(numel(strengthIntraI{IX})+numel(strengthInterI{IX}))*(mean(strengthIntraI{IX})+mean(strengthInterI{IX}))/2/9;
    if ~isnan(WW(IX))
        for k = 1:round(WW(IX))
            L1 = [L1,XX(IX)];
        end
    end
    XX2(IX) = mean(interDiff_exc{IX}) - mean(intraDiff_exc{IX});
    WW2(IX) = sqrt(numel(strengthIntraE{IX})+numel(strengthInterE{IX}))*(mean(strengthIntraE{IX})+mean(strengthInterE{IX}))/2/9;
    if ~isnan(WW2(IX))
        for k = 1:round(WW2(IX))
            L2 = [L2,XX2(IX)];
        end
    end
end

figure(31);  hold on;
% rather than a square plot, make it thinner
violinPlot(L1', 'histOri', 'left', 'widthDiv', [2 1], 'showMM', 0, ...
    'color',  mat2cell(colors(1, : ), 1));
 
violinPlot(L2', 'histOri', 'right', 'widthDiv', [2 2], 'showMM', 0, ...
    'color',  [0 0 0]);
set(gca, 'xtick', [0.6 1.4], 'xticklabel', {'inhibition', 'excitation'}, 'xlim', [0.2 1.8]);
ylabel('tuning strength [pA]'); xlabel('weighted smooth distributions');
 
% add significance stars for each bar
xticks = get(gca, 'xtick');
for b = 1:2,
    if b == 1; XL = L1; else; XL = L2; end
    [~, pval] = ttest(XL);
    yval = max(XL) * 1.2; % plot this on top of the bar
    yval = 6; % plot below
    mysigstar(gca, xticks(b), yval, pval);
    % if mysigstar gets just 1 xpos input, it will only plot stars
end
 
XX = []; YY = []; XX2 = []; YY2 = []; 
for IX = 1:98
    XX(IX) = mean(interDiff_inh{IX}) - mean(intraDiff_inh{IX});
    WW(IX) = sqrt(numel(strengthIntraI{IX})+numel(strengthInterI{IX}))*(mean(strengthIntraI{IX})+mean(strengthInterI{IX}))/2/9;
    XX2(IX) = mean(interDiff_exc{IX}) - mean(intraDiff_exc{IX});
    WW2(IX) = sqrt(numel(strengthIntraE{IX})+numel(strengthInterE{IX}))*(mean(strengthIntraE{IX})+mean(strengthInterE{IX}))/2/9;
%     plot(XX(IX),'.','MarkerSize',max(2,WW(IX))); hold on;
    offset1 = rand()*0.08;
%     if WW(IX) > 10 || WW2(IX) > 10
%     plot([XX(IX),XX2(IX)],[0.25+offset1,0.3+offset1],'Color',[0.8 0.8 0.8],'LineWidth',sqrt(max(0,WW2(IX))*max(0,WW(IX)))/10 ); hold on;
%     end
    plot(0.4+offset1,XX(IX),'.r','MarkerSize',sqrt(max(2,WW(IX))*max(2,WW(IX)))); hold on;
    plot(1.6+offset1,XX2(IX),'.k','MarkerSize',sqrt(max(2,WW2(IX))*max(2,WW2(IX)))); hold on;

%     plot(0:0.1:20,0:0.1:20,'k');
end
hold off




%% select good cells
% IXperm = randperm(98);
for IX = 1:98
%         answer_2{IX} = manualSelection(IX,X_all,datasetSingleCells,days,onsetZ);
end

% post-selection
for IX = 1:98
    if str2num(answer{IX}{1})
        A1 = X_all{IX}.onsetVC70; A1(A1 == 0) = NaN;
        A1 = squeeze(nanmean(A1,2));
        A2 = X_all{IX}.onsetVC00; A2(A2 == 0) = NaN;
        A2 = squeeze( nanmean(A2,2));
        figure(413); plot(A1'); hold on; plot(A2'); hold off;
        answerX = inputdlg('Excellent?','',1);
        answer_2{IX} = str2double(answerX);
    else
        answer_2{IX} = 0;
    end
end

for IX = 1:98
    for k = 2:3
        A = answer{IX}{k};
        D = strsplit(A);
        if ~cellfun(@isempty,D)
            for j = 1:numel(D)
                chosen_ones(IX,k-1,j) = str2num(D{j});
            end
        end
    end
end


