%% go to main folder
cd('C:\Data\rupppete\PhD\electrophysiology2016')

%% show annotation viewer
xx = pwd;
cd('C:\Data\rupppete\PhD\electrophysiology2016\AnnotationViewer');
load('StacksUint16.mat');
figure(99);
imshow3DofDp(Horizontal,Sagittal,'NeuronList_Annotation',[8387 8557]);
cd(xx)

%show Dp stack of this fish
showDpStack()

% show stack of the current cell
showCellStack()

% show recordings for the current cell
showEphys()

