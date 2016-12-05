
%% Dp viewer

load('StacksUint16.mat');

figure(99);

imshow3DofDp(Horizontal,Sagittal,'NeuronList_Annotation',[8387 8557]);


figure(112);

for k = 45:100
    imagesc(Horizontal(:,:,k),[8387 8557]); hold on; colormap(gray)
    for j = 1:numel(NeuronList)
        if abs(NeuronList{j}.z - k) < 2
            plot(NeuronList{j}.x,NeuronList{j}.y,'y','MarkerSize',12);
        end
    end
    hold off;
    pause(0.5)
end

% key short cuts:
x - 
w -
e/q (?) - 
s -
d - 
b -




