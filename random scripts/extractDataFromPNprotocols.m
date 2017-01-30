


 
trialSelect = [13 16 13 16 13 15 13 13 13 13 16 16 19 13 20 14 4  14 4  13 18 16 13 13 14 16 16 17 15 19 14 12  17  19  21 15  10  15  21  13  17];
cellIX = [64 66 67 68 69 70 72 73 74 75 76 77 78 79 80 81 82 83 84 85 86 87 88 89 90 92 93 94 95 96 97 99 100 101 104 107 108 109 111 112 113];
paradigm = [1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  3  3  3  3  2  2  2  2  2  3  3  3  3  3  3  3  1   1   1   1  1   1   1   1  1   1];
for j = 1:numel(trialSelect)
    
    
   
    cellIX(j)
    cd(['C:\Data\rupppete\PhD\electrophysiology2016\SingleCells\',datasetSingleCells{cellIX(j)}.CellID]);
    
    FileList = dir('**.xsg');
    load(FileList(trialSelect(j)).name,'-mat');
    A = data.ephys.trace_1;
    
    if paradigm(j) == 1
        excMean = A(251:150250);
        excMean = reshape(excMean,1500,100);
        
        baseline = excMean(1150:1500,:);
        inline = excMean(200:1000,:);

    elseif paradigm(j) == 2
        excMean = A(10001:210000);
        excMean = reshape(excMean,2000,100);
        
        baseline = excMean(1750:2000,:);
        inline = excMean(500:1500,:);
    else
        excMean = A(10001:199000);
        excMean = reshape(excMean,1890,100);
        
        baseline = excMean(1550:1890,:);
        inline = excMean(500:1400,:);
    end
    
    Rin = -10./( median(baseline) - median(inline) );

    tttrace = median(excMean');

    Rs1 = 10e-3./(median(tttrace)-min(tttrace))*1e12/1e6;
    Rs2 = - 10e-3./(median(tttrace)-max(tttrace))*1e12/1e6;

    Y = median(excMean');
    X = (0:numel(Y)-1)/10;
    diffX = Y - circshift(Y,[0 -1]); startIDX = find(diffX>0,1,'first');
    [xData, yData] = prepareCurveData( X(startIDX:400), Y(startIDX:400) );
    ft = fittype( 'a*exp(-x/b)+d*exp(-x/e)+c', 'independent', 'x', 'dependent', 'y' );
    opts = fitoptions( ft ); opts.Display = 'Off';
    opts.Lower = [-Inf -Inf -Inf -Inf -Inf];
    opts.StartPoint = [115 3 -22 22 1];
    opts.Upper = [Inf Inf Inf Inf Inf];
    try
        [fitresult, gof] = fit( xData, yData, ft, opts );
        tauVC = max(fitresult.b,fitresult.e);
    catch
        tauVC = NaN;
    end
    

    biophysics{cellIX(j)}.Rin = Rin*1e9;
    biophysics{cellIX(j)}.Rs = [Rs1,Rs2];

    biophysics{cellIX(j)}.Cm = tauVC*1e-3/mean([Rs1 Rs2])/1e6;
    biophysics{cellIX(j)}.tau = tauVC*1e-3/mean([Rs1 Rs2])/1e6*median(Rin(:))*1e9;    
    
    
end
    

clear goodX
for j = 1:numel(cellIX)
    [mean(biophysics{cellIX(j)}.Rin(:))/1e9,mean(biophysics{cellIX(j)}.Rs(:)),mean(biophysics{cellIX(j)}.tau(:))*1000,j]
    datasetSingleCells{cellIX(j)}
    goodX{j} = inputdlg({'Good','TimeConstant'},'UI / Cell selection',1);

end

for j = 1:numel(cellIX)
    
    if strcmp(goodX{j}(1),'0')
        biophysics{cellIX(j)}.Rin = NaN;
    end
    if strcmp(goodX{j}(2),'0')
        biophysics{cellIX(j)}.tau = NaN;
        biophysics{cellIX(j)}.Cm = NaN;
    end
end
  

    
    