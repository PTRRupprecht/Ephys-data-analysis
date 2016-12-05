
% get data points that define the surface of the brain (manual annotation
% using the annotation viewer)
clear x y z 
for j =1:numel(NeuronList)% 1:259
    x(j) = NeuronList{j}.x;
    y(j) = NeuronList{j}.y;
    z(j) = NeuronList{j}.z;
end

% save('plane2.mat','x','y','z')


%% Fit a plane to the datapoints (parameter choice done using cftool)
[xData, yData, zData] = prepareSurfaceData( x, y, z );
ft = fittype( 'loess' ); opts = fitoptions( ft ); opts.Robust = 'Bisquare';
opts.Span = 0.05; opts.Normalize = 'on';
[fitresultLoess, gof] = fit( [xData, yData], zData, ft, opts );

%% generate projection of zstacks according to the boundary plane
z_coord = zeros(1490,787);
[Xgrid, Ygrid] = meshgrid(1:1490,1:787);
Zgrid = fitresultLoess(Xgrid,Ygrid);

Map = zeros(size(Zgrid)); % horizontal stack
% MapY = zeros(size(Zgrid)); % sagittal stack
for j = 1:size(Zgrid,1)
    for k = 1:size(Zgrid,2)
        Map(j,k) = Horizontal(j,k,round(min(max(Zgrid(j,k),1),291)));
%         MapY(j,k) = Sagittal(round(min(max(1,Zgrid(j,k)*6),1746)),k,ceil(j/6));
    end
end
% show projections
figure(3); imagesc(Map,[8400 8590]); colormap(gray)
hold on;
cmap = jet(numel(xyz(:,1)));
for j = 1:size(xyz,1)
    plot(xyz(j,1),xyz(j,2),'.','MarkerSize',28,'Color',cmap(j,:)); hold on;
end


%% part II of the projection : sagittal surface projection
opts.Span = 0.05; opts.Normalize = 'on';
[fitresultLoess, gof] = fit( [xData, zData*6], yData, ft, opts );

%% generate projection of zstacks according to the boundary plane
z_coord = zeros(1490,1080);
[Xgrid, Ygrid] = meshgrid(1:1490,1:1080);
Zgrid = fitresultLoess(Xgrid,Ygrid);

Zgrid(Zgrid>1080) = 1080;
Zgrid(Zgrid<1) = 1;

MapY = zeros(size(Zgrid)); % sagittal stack
for j = 1:size(Zgrid,1)
    for k = 1:size(Zgrid,2)
        MapY(j,k) = Sagittal(j,k,round(min(max(1,Zgrid(j,k)/6),1080)));
    end
end
% show projection
figure(2); imagesc(MapY,[8400 8590]); colormap(gray)
hold on;
cmap = jet(numel(xyz(:,1)));
for j = 1:size(xyz,1)
    plot(xyz(j,1),xyz(j,3)*6,'.','MarkerSize',28,'Color',cmap(j,:)); hold on;
end




%% correct for stretching in y (approximately)

% 'newpixels' gives a mapping of input pixels in y (1:787) to the target
% pixels in y (1:2400); the factor 6 accounts for the fact that resolution
% in z is 6x lower than in x/y

testY = Zgrid(:,1200); % this is a good proxy for overall stretching
newGrid = cumsum(1+6*abs(smooth(diff(testY,1),50)));
newGrid = newGrid - min(newGrid) + 1;
newpixels = interp1(newGrid,1:786,1:round(max(newGrid)));

% generate the grid points along the boundary plane; these define
[Xgrid, Ygrid] = meshgrid(1:1490,780:length(newpixels)); % 780 = arbitrary, after that, there is no brain tissue
Zgrid = fitresultLoess(Xgrid,newpixels(Ygrid));

% generate the (streched) projections of the boundary plane to 2D
MapX = zeros(size(Zgrid));
% MapX2 = zeros(size(Zgrid));
MapY = zeros(size(Zgrid));
for j = 1:size(Zgrid,1)
    if mod(j,50) ==0; disp(j/size(Zgrid,1)); end
    for k = 1:size(Zgrid,2)
        newj = round(interp1(1:numel(newpixels),newpixels,j+779));
        ixZ = min(max(Zgrid(j,k),1),291);
        ixZ2 = min(max(1,Zgrid(j,k)*6),1746);
        MapX(j,k) = Horizontal(newj,k,round(ixZ)) ;
%         MapX2(j,k) = Horizontal(newj,k,floor(ixZ))*(ixZ-floor(ixZ)) + Horizontal(newj,k,ceil(ixZ))*(ceil(ixZ)-ixZ) ;
        MapY(j,k) = Sagittal(round(ixZ2),k,ceil(newj/6));
    end
end
figure(23); imagesc(MapX,[8400 8590]); colormap(gray)
figure(22); imagesc(MapY,[8400 8590]); colormap(gray)
% combine the horizontal and the sagittal stack based on a manually
% selected mask (BW)

% BW = roipoly;
BWsmooth = conv2(double(BW),fspecial('gaussian',[25 25], 25),'same');
Map_combined = MapX.*(1-BWsmooth) + BWsmooth.*MapY;

figure(41), imagesc(Map_combined,[8410 8570]); colormap(gray)

%% contour overlay

ZgridO = Zgrid; ZgridO(ZgridO>180) = 180;
hold on; contour(ZgridO,[0:10:140],'Color',[1.0 0.7 0.7])


%% load point cloud from patching dataset

datasetList = dir('dataset*.mat');
load(datasetList(1).name);
xyz = zeros(numel(datasetSingleCells),3);
for k = 1:numel(datasetSingleCells)
    xyz(k,:) = datasetSingleCells{k}.pos;
end


%% find closest corresponding point on surface, given a point

for j = 1:size(xyz,1)
    if ~isnan(xyz(j,1))
        pp = xyz(j,:);
        pp(2) = interp1(newpixels(780:length(newpixels)),(780:length(newpixels))-779,pp(2)); % account for distortions in y
        pp(3) = pp(3)*6;
        [px,py] = meshgrid(1:size(ZgridO,1),1:size(ZgridO,2));
        pz = ZgridO'*6;

        distanceGrid = zeros(size(py,1),size(py,2));
        for jj = 1:size(py,1);
            ww = [py(jj,1:size(py,2)); px(jj,1:size(py,2)); pz(jj,1:size(py,2))];
            distanceGrid(jj,1:size(py,2)) = sum((ww- repmat(pp',[1 size(py,2)])).^2);
        end

        minX = min(distanceGrid(:));
        [gx(j),gy(j)] = find(distanceGrid == minX);
        hold on; plot(gx(j),gy(j),'.k','Color',cmap(j,:),'MarkerSize',28)
        pause(0.1)
    end
end


%% manually select neurons using the map
[xx,yy] = ginput(1);
distanceVector = sum(([gx;gy]- repmat([xx;yy],[1 size(gy,2)])).^2);
[~,IX] = min(distanceVector);

datasetSingleCells{IX}



