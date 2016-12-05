%% temporal delay between the three odor lines

% Dp neurons odor onset

Dp = [ 12.09 12.3 12.59;
    11.85 11.9 12.3;
    11.72 11.67 12.1];
X = Dp;
delay1 = X(1,:) - X(3,:);
delay2 = X(2,:) - X(3,:);

OB = [ 12.31 12.13 12.75 12.36 12.25 12.26 12.32 12.28 12.36;
    11.96 11.9 NaN 11.95 NaN NaN NaN 11.91 11.96;
    11.67 11.56 11.88 11.69 11.66 11.69 11.62 11.52 11.67];

X = OB;
delay1o = X(1,:) - X(3,:);
delay2o = X(2,:) - X(3,:);

[ nanmean([delay1o, delay1]) nanstd([delay1o, delay1])  nanstd([delay1o, delay1])/sum(~isnan([delay1o, delay1]),2) ]
[ nanmean([delay2o, delay2]) nanstd([delay2o, delay2])  nanstd([delay2o, delay2])/sum(~isnan([delay2o, delay2]),2) ]


% result: mean +- std error (only based on OB data)

delay1 = (673 +- 11) ms  (variability of measurement std = 98 ms)
delay2 = (314 +- 10) ms  (variability of measurement std = 51 ms)


save(['delays',date,'.mat'],'delay1','delay2','delay3')


% update the dataset

index = 98
datasetSingleCells{index}.odors = {'Food' 'Arg' 'TDC'};
datasetSingleCells{index}.odorLine = [1 2 3];
index+1
datasetSingleCells{index+1}


datasetSingleCells{index}



