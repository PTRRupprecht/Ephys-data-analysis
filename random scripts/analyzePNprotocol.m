
%% open all P/N protocols

clear all
load('dataset27-Jul-2016_80.mat');

counter = 1;
for k = 1:numel(datasetSingleCells)
    try
        if isempty(datasetSingleCells{k}.PN70) == 0
            numel(datasetSingleCells{k}.PN70)
            for jjj = 1:numel(datasetSingleCells{k}.PN70)
                index = datasetSingleCells{k}.PN70(jjj);
                cd(datasetSingleCells{k}.CellID)
                listXSG = dir('*.xsg');
                load(listXSG(index).name,'-mat');
                A = data.ephys.trace_1;

                transients = zeros(1500,400);
                for kk = 1:400
                    transients(:,kk) = A((1:1500)+(kk-1)*1500);
                end

                figure(40);
                timecourse2D = mean(transients(:,1:100),2);
                subplot(4,5,counter); hold on;
                plot((1:1500)/10,timecourse2D);
                axis([22 45 -1500 500]);
                xlabel(strcat(datasetSingleCells{k}.CellID,32,32,num2str(round(10000/(max(timecourse2D) - mean(timecourse2D(1:25)))))),'Interpreter','None')
                 
                
                sum(timecourse2D(1:1000) - mean(timecourse2D(1:250)))
                
                
                figure(41);
                timecourse2D = transients(:,101:200) - repmat(mean(transients(:,1:100),2)*7,[1 100]);
                subplot(4,5,counter); hold on; cmap = jet(100);
                for kk = 1:100
                    plot((1:1500)/10,smooth(timecourse2D(:,kk),1),'Color',cmap(kk,:));
                end
                axis([22 45 -1500 500]);
                xlabel(datasetSingleCells{k}.CellID,'Interpreter','None')

                figure(42);
                timecourse2D = transients(:,301:400) - repmat(mean(transients(:,201:300),2)*3.5,[1 100]);
                subplot(4,5,counter); hold on; cmap = jet(100);
                for kk = 1:100
                    plot((1:1500)/10,smooth(timecourse2D(:,kk),1),'Color',cmap(kk,:));
                end
                axis([22 45 -1500 500]);
                xlabel(datasetSingleCells{k}.CellID,'Interpreter','None')
                counter = counter + 1;
                drawnow;
                cd ..
            end
        end
    end
end

