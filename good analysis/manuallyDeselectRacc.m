
function stableResponse = manuallyDeselectRacc(IX,datasetSingleCell,sTrace,ssTrace)
%manualSelection(IX,X_all,datasetSingleCells,days,onsetZ)

figure(3); 
offsetVC = 0; cmap = distinguishable_colors(3);
for k = 1:size(sTrace.onsetVC70,1)
    subplot(1,2,1);
    for jj = 1:size(sTrace.onsetVC70,2)
        AA = squeeze(sTrace.onsetVC70(k,jj,1:20:end));
        plot((1:20:6e5)/1e4,AA(1:1:end)+offsetVC,'Color',cmap(k,:)); hold on;
        text(4.1e5/1e4,nanmedian(AA)+offsetVC,num2str(ssTrace.trialsVC70(k,jj)),'FontSize',12)
        if ~isnan(AA(1))
            offsetVC = offsetVC - min(AA);
        end
    end
    xlim([0 38]);
end
hold off; 
title('Excitatory currents');


figure(3); offsetVC = 0;
for k = 1:size(sTrace.onsetVC00,1)
    k
    subplot(1,2,2);
    for jj = 1:size(sTrace.onsetVC00,2)
        AA = (squeeze(sTrace.onsetVC00(k,jj,1:5:end)));
        plot((1:20:6e5)/1e4,AA(1:4:end)+offsetVC,'Color',cmap(k,:)); hold on;
        text(4.1e5/1e4,nanmedian(AA)+offsetVC,num2str(ssTrace.trialsVC00(k,jj)),'FontSize',12)
        if ~isnan(AA(1))
            offsetVC = offsetVC + max(AA);
        end
    end
    xlim([0 38]);
end
hold off; 
title('Inhibitory currents');

stableResponse = inputdlg({'Enter good trials, exc:','Enter good trials, inh:'},'UI for Trial selection',1);


%% load raw file and calculate cell parameters from transients
  

end



