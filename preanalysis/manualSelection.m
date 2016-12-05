
function answer = manualSelection(IX,X_all,datasetSingleCells,days,onsetZ)

    % IX = 76;
    AA = X_all{IX}.onsetVC00;
    
    figure(22); subplot(2,1,1); hold off; plot(0,0); subplot(2,1,2);hold off; plot(0,0);
    cmap = distinguishable_colors(13);
    for k = 1:size(AA,1)
        figure(2); subplot(1,3,k); hold off; plot(0,0);
        figure(2); subplot(1,3,k);
        plot([10 10],[-50 1300],'k'); hold on;
        plot([12 12],[-50 1300],'k');
        offsetV = 0;
        odorIndex = find(~cellfun(@isempty,strfind(datasetSingleCells{IX}.odors,X_all{IX}.odorsVC{k} ) ));
        trials_T = find(~cellfun(@isempty,strfind(datasetSingleCells{IX}.VC0odor,X_all{IX}.odorsVC{k}) ));
        trials = datasetSingleCells{IX}.VC0(trials_T);
        dayIndex = find(~cellfun(@isempty,strfind(days,X_all{IX}.dateID)));
        delay = onsetZ(dayIndex,odorIndex);
        for j = 1:numel(trials)
            D = squeeze(circshift(AA(k,j,:),[0 0 1e5-round(delay)]));
            if ~isnan(sum(D))
                L_trace = D;
                xdata = [1:size(L_trace)]/1e4; %in sec
                new_data = timeseries(L_trace,xdata);
                idealfilter_data = idealfilter(new_data,[5 200],'pass');
                std_noise = std(idealfilter_data.getdatasamples(1:0.95e4));
                std_top = quantile(idealfilter_data.getdatasamples(1:0.95e4),0.02);
                std_low = quantile(idealfilter_data.getdatasamples(1:0.95e4),0.98);
                std_noise = std_top - std_low;
                AUC = nanmean(squeeze(D(1e5-500:1.2e5))) - nanmean(squeeze(D(1:0.95e5)));
                figure(2); subplot(1,3,k); hold on; plot( (1:numel(D))/1e4,smooth(D,100)+offsetV,'Color',cmap(k,:));
                text(numel(D)/1e4+1,median(D(1:10:end)+offsetV),[num2str(trials(j)),',',32,num2str(round(AUC)),',',32,num2str(round(std_noise))],'FontSize',12)
                offsetV = offsetV + max(smooth(D,50))-min(smooth(D,50));
                ylim([-50 offsetV+100]);
                figure(22); subplot(2,1,1); hold on; plot(trials(j),AUC,'.','Color',cmap(k,:),'MarkerSize',17);
                subplot(2,1,2); hold on; plot(trials(j),std_noise,'.','Color',cmap(k,:),'MarkerSize',17);
            end
        end
        %     axis off
    end
    for k = size( AA,1)+1:3
        figure(2); subplot(1,3,k);  hold off; plot(1,1);
    end
    
    
    AA = X_all{IX}.onsetVC70;
    
    cmap = distinguishable_colors(13);
    for k = 1:size(AA,1)
        figure(3); subplot(1,3,k); hold off; plot(0,0);
        figure(3);subplot(1,3,k);
        plot([10 10],[-150 1300],'k'); hold on;
        plot([12 12],[-150 1300],'k');offsetV = 0;
        odorIndex = find(~cellfun(@isempty,strfind(datasetSingleCells{IX}.odors,X_all{IX}.odorsVC{k} ) ));
        trials_T = find(~cellfun(@isempty,strfind(datasetSingleCells{IX}.VC70odor,X_all{IX}.odorsVC{k}) ));
        trials = datasetSingleCells{IX}.VC70(trials_T);
        dayIndex = find(~cellfun(@isempty,strfind(days,X_all{IX}.dateID)));
        delay = onsetZ(dayIndex,odorIndex);
        for j = 1:numel(trials)
            D = squeeze(circshift(AA(k,j,:),[0 0 1e5-round(delay)]));
            if ~isnan(sum(D))
                L_trace = D;
                xdata = [1:size(L_trace)]/1e4; %in sec
                new_data = timeseries(L_trace,xdata);
                idealfilter_data = idealfilter(new_data,[5 200],'pass');
                std_noise = std(idealfilter_data.getdatasamples(1:0.95e4));
                std_top = quantile(idealfilter_data.getdatasamples(1:0.95e4),0.98);
                std_low = quantile(idealfilter_data.getdatasamples(1:0.95e4),0.02);
                std_noise = std_top - std_low;
                AUC = nanmean(squeeze(D(1e5-500:1.2e5))) - nanmean(squeeze(D(1:0.95e5)));
                figure(3); subplot(1,3,k); hold on; plot( (1:numel(D))/1e4,smooth(D,100)+offsetV,'Color',cmap(k,:)); hold on;
                text(numel(D)/1e4+1,median(D(1:10:end)+offsetV),[num2str(trials(j)),',',32,num2str(round(AUC)),',',32,num2str(round(std_noise))],'FontSize',12)
                offsetV = offsetV + max(smooth(D,50))-min(smooth(D,50));
                ylim([-150 offsetV+100]);
                figure(22); subplot(2,1,1); hold on; plot(trials(j),AUC,'.','Color',cmap(k,:),'Marker','x','MarkerSize',8,'LineWidth',2);
                subplot(2,1,2); hold on; plot(trials(j),std_noise,'.','Color',cmap(k,:),'Marker','x','MarkerSize',8,'LineWidth',2);
            end
        end
        %     axis off
    end
    if size( AA,1)<3
        for k = size( AA,1)+1:3
        figure(3); subplot(1,3,k); hold off; plot(1,1);
        end
    end
    
    datasetSingleCells{IX}
    
    NeuronViewerMap(1,IX,4831);

    answer = inputdlg({'Transients:','Enter good trials, inh:','Enter good trials, exc:'},'User Interface for Trial selection',1);

end



