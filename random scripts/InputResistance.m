
counter = 1;
for k = [93 95 96 97 92 86 83 82 81]
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
        [xData, yData] = prepareCurveData( [],  mean(transients{1}(103:250,:),2) );
        ft = fittype( 'a*exp(-x/b)+c', 'independent', 'x', 'dependent', 'y' );
        opts = fitoptions( ft );
        opts.Display = 'Off';
        opts.Lower = [-Inf -Inf -Inf];
        opts.StartPoint = [20 40 0.48732168329619];
        opts.Upper = [Inf Inf Inf];
        [fitresult, gof] = fit( xData, yData, ft, opts );
        timeConst(counter) = fitresult.b/10';
        
        
        R1(counter) = (mean(mean(transients{1}(1200:1400,:),2),1) -  mean(mean(transients{1}(1:90,:),2),1))/10;
        R2(counter) = (mean(mean(transients{2}(1200:1400,:),2),1) -  mean(mean(transients{2}(1:90,:),2),1))/20;
        R3(counter) = (mean(mean(transients{3}(1200:1400,:),2),1) -  mean(mean(transients{3}(1:90,:),2),1))/30;

        Cm(counter) = timeConst(counter)/R1(counter)/1e6;
        
%         R1(counter) = (mean(mean(transients{1}(1200:1400,:),2),1) -  mean(A(20.5e4:24e4)))/10;
%         R2(counter) = (mean(mean(transients{2}(1200:1400,:),2),1) -  mean(A(20.5e4:24e4)))/20;
%         R3(counter) = (mean(mean(transients{3}(1200:1400,:),2),1) -  mean(A(20.5e4:24e4)))/30;
counter = counter + 1;
    end
    cd ..
end


