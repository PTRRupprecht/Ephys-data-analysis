
%% open all P/N protocols

clear all
load('dataset01-Aug-2016_97.mat');

counter = 1;
for k = 1:98%numel(datasetSingleCells)
    try
        if isempty(datasetSingleCells{k}.PN70II) == 0
            numel(datasetSingleCells{k}.PN70II)
            for jjj = 1:numel(datasetSingleCells{k}.PN70II)
                index = datasetSingleCells{k}.PN70II(jjj);
                cd(datasetSingleCells{k}.CellID)
                listXSG = dir('*.xsg');
                load(listXSG(index).name,'-mat');
                A = data.ephys.trace_1;

                if k < 90 && k > 84
                    duration = 2000;
                else
                    duration = 1890;
                end
                
                transients{1} = zeros(duration,100);
                for kk = 1:100
                    transients{1}(:,kk) = A((1:duration)+(kk-1)*duration + 1e4-100);
                end
                for j = 2:7
                    transients{j} = zeros(duration,30);
                    for kk = 1:30
                        transients{j}(:,kk) = A((1:duration)+(kk-1)*duration+1e5*j+5e4-100);
                    end
                end
                figure(32+k*10+jjj); hold on;
                for j = 2:7
                    subplot(2,3,j-1); hold on;
                    cmap = paruly(30);
                    temp = zeros(duration,1);
                    for kj = 1:30
                        plot((1:duration)/10,smooth(transients{j}(:,kj)'-j*mean(transients{1}'),8),'Color',cmap(kj,:))
                        temp = temp + transients{j}(:,kj)-j*mean(transients{1}')';
                        axis([0 duration/10 -500 300])
                    end
                    plot((1:duration)/10,temp/30,'Color',[0 0 0])
                    xlabel('time [ms]'); ylabel('current, P/N-corrected [pA]')
                end
                suptitle([datasetSingleCells{k}.CellID,32,'trial',32,num2str(jjj),32,'cell no.',32,num2str(k)])
            end
        end
        cd ..
    catch
        cd('C:\Data\rupppete\PhD\electrophysiology2016\SingleCells');
        disp(['Fail!',32,num2str(k)])
    end
end

