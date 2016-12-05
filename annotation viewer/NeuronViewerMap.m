function NeuronViewerMap(value,neuronID,figureNb)

colors = numel(unique(value));
if colors > 10
    cmap = brewermap(256,'OrRd');
    value = round((value - min(value))/(max(value)- min(value))*255) + 1;
else
    cmap = distinguishable_colors(colors);%     brewermap(colors,'OrRd');% winter(colors);
    if min(value) == 0; 
        value = value + 1;
    end
end

tempFolder = pwd;

cd('C:\Data\rupppete\PhD\electrophysiology2016\SingleCells');
datasetList = dir('dataset*.mat');
load(datasetList(1).name);

xyz = zeros(numel(datasetSingleCells),3);
for k = 1:numel(datasetSingleCells)
    xyz(k,:) = datasetSingleCells{k}.pos;
end


cd('C:\Data\rupppete\PhD\electrophysiology2016\AnnotationViewer');
load('BoundaryMapping.mat');
load('SimpleMapsXY.mat');
load('coordinates_07-Nov-2016.mat');

Map(Map>8700) = 8700;
figure(figureNb+1); imagesc(Map,[8380 8790]); colormap(gray); axis equal off
hold on;
for k = 1:numel(value)
    try
        plot(xyz(neuronID(k),1)+round((rand(1)-0.5)*30),xyz(neuronID(k),2)+round((rand(1)-0.5)*30),'.','MarkerSize',25,'Color',cmap(round(value(k)),:));
    catch
        k
        round(value(k))
        plot(100+round((rand(1)-0.5)*30),100+round((rand(1)-0.5)*30),'.','MarkerSize',5,'Color',cmap(round(value(k)),:));
    end
end
MapY(MapY>8700) = 8700;
figure(figureNb+2); imagesc(MapY,[8380 8790]); colormap(gray); axis equal off
hold on;
for k = 1:numel(value)
    try
        plot(xyz(neuronID(k),1)+round((rand(1)-0.5)*30),xyz(neuronID(k),3)*6+round((rand(1)-0.5)*30),'.','MarkerSize',25,'Color',cmap(round(value(k)),:));
    catch
        k
        round(value(k))
        plot(100+round((rand(1)-0.5)*30),100+round((rand(1)-0.5)*30),'.','MarkerSize',5,'Color',cmap(round(value(k)),:));
    end
end

Map_combined(Map_combined>8700) = 8700;
figure(figureNb), imagesc(Map_combined,[8380 8790]); colormap(gray); axis equal off;
% akZoom_PR()
ZgridO = Zgrid; ZgridO(ZgridO>180) = 180;
hold on; contour(ZgridO,[0:10:140],'Color',[1 0.8 0.9])
rng(7831)
for k = 1:numel(value)
    try
        plot(gx(neuronID(k))+round((rand(1)-0.5)*30),gy(neuronID(k))+round((rand(1)-0.5)*30),'.','MarkerSize',25,'Color',cmap(round(value(k)),:));
    catch
        k
        round(value(k))
        plot(100+round((rand(1)-0.5)*30),100+round((rand(1)-0.5)*30),'.','MarkerSize',25,'Color',cmap(round(value(k)),:));
    end
end


% figure(44), imagesc(Map_combined,[8410 8570]); colormap(gray); axis equal off;
% akZoom_PR()
% ZgridO = Zgrid; ZgridO(ZgridO>180) = 180;
% hold on; contour(ZgridO,[0:10:140],'Color',[1 0.8 0.9])
% for k = 1:numel(gx)
%     if strcmp(datasetSingleCells{k}.CellID(1:6),'160731')
%         plot(gx(k),gy(k),'.','MarkerSize',20,'Color',cmap(values(k),:));
%     end
% end
% set(gcf, 'WindowKeyPressFcn', {@showRawData,gx,gy});
    
cd(tempFolder)
