function  imshow3DofDp(Horizontal,Sagittal,outputfile,disprange)
%imshow3Dfull1 displays 3D grayscale images from three perpendicular views
%(i.e. axial, sagittal, and coronal) in slice by slice fashion with mouse
%based slice browsing and window and level adjustment control.
%
% Usage:
% imshow3Dfull ( Image )
% imshow3Dfull ( Image , [] )
% imshow3Dfull ( Image , [LOW HIGH] )
%   
%    Image:      3D image MxNxK (K slices of MxN images) 
%    [LOW HIGH]: display range that controls the display intensity range of
%                a grayscale image (default: the widest available range)
%
% Use the scroll bar or mouse scroll wheel to switch between slices. To
% adjust window and level values keep the mouse right button pressed and
% drag the mouse up and down (for level adjustment) or right and left (for
% window adjustment).
%
% Use 'A', 'S', and 'C' buttons to switch between axial, sagittal and
% coronal views, respectivelly.
% 
% "Auto W/L" button adjust the window and level automatically 
%
% While "Fine Tune" check box is checked the window/level adjustment gets
% 16 times less sensitive to mouse movement, to make it easier to control
% display intensity rang.
%
% Note: The sensitivity of mouse based window and level adjustment is set
% based on the user defined display intensity range; the wider the range
% the more sensitivity to mouse drag.
% 
% 
%   Example
%   --------
%       % Display an image (MRI example)
%       load mri 
%       Image = squeeze(D); 
%       figure, 
%       imshow3Dfull(Image) 
%
%       % Display the image, adjust the display range
%       figure,
%       imshow3Dfull(Image,[20 100]);
%
%   See also IMSHOW.

%
% - Maysam Shahedi (mshahedi@gmail.com)
% - Released: 1.0.0   Date: 2013/04/15
% - Revision: 1.1.0   Date: 2013/04/19
% - Revision: 2.0.0   Date: 2014/08/05
% - Peter Rupprecht (p.t.r.rupprecht@gmail.com)
% - Revision: 1.0     Date: 2016/04/04

sno = size(Horizontal);  % image size
sno_a = sno(3);  % number of axial slices
S_a = round(sno_a/2);
sno_s = sno(2);  % number of sagittal slices
S_s = round(sno_s/2);
sno_c = sno(1);  % number of coronal slices
S_c = round(sno_c/2);
S = S_a;
sno = sno_a;


global InitialCoord;
global NeuronList;
global blobs;
blobs = 0;
viewport = 1;
NeuronList(:) = [];

MinV = 0;
MaxV = max(Horizontal(:));
LevV = (double( MaxV) + double(MinV)) / 2;
Win = double(MaxV) - double(MinV);
WLAdjCoe = (Win + 1)/1024;
FineTuneC = [1 1/16];    % Regular/Fine-tune mode coefficients

if isa(Horizontal,'uint8')
    MaxV = uint8(Inf);
    MinV = uint8(-Inf);
    LevV = (double( MaxV) + double(MinV)) / 2;
    Win = double(MaxV) - double(MinV);
    WLAdjCoe = (Win + 1)/1024;
elseif isa(Horizontal,'uint16')
    MaxV = uint16(Inf);
    MinV = uint16(-Inf);
    LevV = (double( MaxV) + double(MinV)) / 2;
    Win = double(MaxV) - double(MinV);
    WLAdjCoe = (Win + 1)/1024;
end 

% ImgAx = Horizontal;
% ImgSg = (permute(Horizontal, [3 1 2]));   % Sagittal view image
% ImgCr = (permute(Horizontal, [3 2 1]));   % Coronal view image

View = 'A';

SFntSz = 9;
LFntSz = 10;
WFntSz = 10;
VwFntSz = 10;
LVFntSz = 9;
WVFntSz = 9;
BtnSz = 10;
ChBxSz = 10;

% disprange = [8496 8496+161];
if (nargin < 2)
    Rmin = 8496;
    Rmax = 161;
%     [Rmin Rmax] = WL2R(Win, LevV);
elseif numel(disprange) == 0
    [Rmin Rmax] = WL2R(Win, LevV);
else
    LevV = (double(disprange(2)) + double(disprange(1))) / 2;
    Win = double(disprange(2)) - double(disprange(1));
    WLAdjCoe = (Win + 1)/1024;
    [Rmin Rmax] = WL2R(Win, LevV);
end

hdl_im = axes('position',[0,0.2,1,0.8]);
imshow(Horizontal(:,:,S), [Rmin Rmax])

FigPos = get(gcf,'Position');
S_Pos = [50 70 uint16(FigPos(3)-100)+1 20];
Stxt_Pos = [50 90 uint16(FigPos(3)-100)+1 15];
Wtxt_Pos = [20 20 60 20];
Wval_Pos = [75 20 60 20];
Ltxt_Pos = [140 20 45 20];
Lval_Pos = [180 20 60 20];
BtnStPnt = uint16(FigPos(3)-210)+1;
if BtnStPnt < 880
    BtnStPnt = 880;
end
Btn_Pos = [BtnStPnt 20 80 20];
ChBx_Pos = [BtnStPnt+90 20 100 20];
Vwtxt_Pos = [255 20 35 20];
VAxBtn_Pos = [290 20 15 20];
VSgBtn_Pos = [310 20 15 20];
VCrBtn_Pos = [330 20 15 20];

VC1 = [380 20 45 20];
VC2 = [430 20 45 20];
VC3 = [480 20 80 20];
VC4 = [565 20 80 20];
VCX = [650 20 500 20];

if sno > 1
    shand = uicontrol('Style', 'slider','Min',1,'Max',sno,'Value',S,'SliderStep',[1/(sno-1) 10/(sno-1)],'Position', S_Pos,'Callback', {@SliceSlider, Horizontal});
    stxthand = uicontrol('Style', 'text','Position', Stxt_Pos,'String',sprintf('Slice# %d / %d',S, sno), 'BackgroundColor', [0.8 0.8 0.8], 'FontSize', SFntSz);
else
    stxthand = uicontrol('Style', 'text','Position', Stxt_Pos,'String','2D image', 'BackgroundColor', [0.8 0.8 0.8], 'FontSize', SFntSz);
end    
ltxthand = uicontrol('Style', 'text','Position', Ltxt_Pos,'String','Level: ', 'BackgroundColor', [0.8 0.8 0.8], 'FontSize', LFntSz);
wtxthand = uicontrol('Style', 'text','Position', Wtxt_Pos,'String','Window: ', 'BackgroundColor', [0.8 0.8 0.8], 'FontSize', WFntSz);
lvalhand = uicontrol('Style', 'edit','Position', Lval_Pos,'String',sprintf('%6.0f',LevV), 'BackgroundColor', [1 1 1], 'FontSize', LVFntSz,'Callback', @WinLevChanged);
wvalhand = uicontrol('Style', 'edit','Position', Wval_Pos,'String',sprintf('%6.0f',Win), 'BackgroundColor', [1 1 1], 'FontSize', WVFntSz,'Callback', @WinLevChanged);
Btnhand = uicontrol('Style', 'pushbutton','Position', Btn_Pos,'String','Auto W/L', 'FontSize', BtnSz, 'Callback' , @AutoAdjust);
ChBxhand = uicontrol('Style', 'checkbox','Position', ChBx_Pos,'String','Fine Tune', 'BackgroundColor', [0.8 0.8 0.8], 'FontSize', ChBxSz);
Vwtxthand = uicontrol('Style', 'text','Position', Vwtxt_Pos,'String','View: ', 'BackgroundColor', [0.8 0.8 0.8], 'FontSize', LFntSz);
VAxBtnhand = uicontrol('Style', 'pushbutton','Position', VAxBtn_Pos,'String','H', 'FontSize', BtnSz, 'Callback' , @AxialView);
VSgBtnhand = uicontrol('Style', 'pushbutton','Position', VSgBtn_Pos,'String','C', 'FontSize', BtnSz, 'Callback' , @SagittalView);
VCrBtnhand = uicontrol('Style', 'pushbutton','Position', VCrBtn_Pos,'String','S', 'FontSize', BtnSz, 'Callback' , @CoronalView);

VCrBtnhand1 = uicontrol('Style', 'pushbutton','Position', VC1,'String','Load', 'FontSize', BtnSz, 'Callback' , @loadX);
VCrBtnhand2 = uicontrol('Style', 'pushbutton','Position', VC2,'String','Save', 'FontSize', BtnSz, 'Callback' , @saveX);
VCrBtnhand3 = uicontrol('Style', 'pushbutton','Position', VC3,'String','Create new', 'FontSize', BtnSz, 'Callback' , @createX);
VCrBtnhand4 = uicontrol('Style', 'pushbutton','Position', VC4,'String','Delete last', 'FontSize', BtnSz, 'Callback' , @deleteX);
Commhand = uicontrol('Style', 'edit','Position', VCX,'String','No comment', 'BackgroundColor', [1 1 1], 'FontSize', LVFntSz);


set (gcf, 'WindowScrollWheelFcn', @mouseScroll);
set (gcf, 'ButtonDownFcn', @mouseClick);
set(get(gca,'Children'),'ButtonDownFcn', @mouseClick);
set(gcf,'WindowButtonUpFcn', @mouseRelease)
set(gcf,'ResizeFcn', @figureResized)

set(gcf, 'WindowKeyPressFcn', {@keypresser});


% -=< Figure resize callback function >=-
    function figureResized(object, eventdata)
        FigPos = get(gcf,'Position');
        S_Pos = [50 45 uint16(FigPos(3)-100)+1 20];
        Stxt_Pos = [50 65 uint16(FigPos(3)-100)+1 15];
        BtnStPnt = uint16(FigPos(3)-210)+1;
        if BtnStPnt < 360
            BtnStPnt = 360;
        end
        Btn_Pos = [BtnStPnt 20 80 20];
        ChBx_Pos = [BtnStPnt+90 20 100 20];
        if sno > 1
            set(shand,'Position', S_Pos);
        end
        set(stxthand,'Position', Stxt_Pos);
        set(ltxthand,'Position', Ltxt_Pos);
        set(wtxthand,'Position', Wtxt_Pos);
        set(lvalhand,'Position', Lval_Pos);
        set(wvalhand,'Position', Wval_Pos);
        set(Btnhand,'Position', Btn_Pos);
        set(ChBxhand,'Position', ChBx_Pos);
        set(Vwtxthand,'Position', Vwtxt_Pos);
        set(VAxBtnhand,'Position', VAxBtn_Pos);
        set(VSgBtnhand,'Position', VSgBtn_Pos);
        set(VCrBtnhand,'Position', VCrBtn_Pos);
        set(VCrBtnhand1,'Position', VC1);
        set(VCrBtnhand2,'Position', VC2);
        set(VCrBtnhand3,'Position', VC3);
        set(VCrBtnhand4,'Position', VC4);
        set(Commhand,'Position', VCX);
    end

% -=< Slice slider callback function >=-
    function SliceSlider (hObj,event, Horizontal)
        S = round(get(hObj,'Value'));
        
        
        if strcmp(View,'A')
            if viewport == 1
                slice = squeeze(Horizontal(:,:,S));
            else
                slice = squeeze(Sagittal(round(S*6),:,:))';
            end
        elseif strcmp(View,'S')
            if viewport == 1
                slice = squeeze(Horizontal(:,S,:))';
            else
                slice = squeeze(Sagittal(:,S,:));
            end
        else
            if viewport == 1
                slice = squeeze(Horizontal(S,:,:));
            else
                slice = squeeze(Sagittal(:,:,round(S/6)));
            end
        end
        if size(slice,1) < 400
            temp = imresize(slice, [size(slice,1)*6.0 size(slice,2)]);
        elseif size(slice,2) < 400
            temp = imresize(slice, [size(slice,1) size(slice,2)*6.0]);
        else
            temp = slice;
        end
            
        K = get(gca,'children');
        set(K(end),'cdata',temp);
        delete(K(1:(end-1)));
        caxis([Rmin Rmax])
        if sno > 1
            set(stxthand, 'String', sprintf('Slice# %d / %d',S, sno));
        else
            set(stxthand, 'String', '2D image');
        end
        if blobs == 1
            cmap = lines(2560); hold on;
            for k = 1:numel(NeuronList)
                if strcmp(View,'A')
                    coox = NeuronList{k}.x; cooy = NeuronList{k}.y; cooz = NeuronList{k}.z; 
                elseif strcmp(View,'S')
                    coox = NeuronList{k}.y; cooy = NeuronList{k}.z*6; cooz = NeuronList{k}.x; 
                else
                    coox = NeuronList{k}.x; cooy = NeuronList{k}.z*6; cooz = NeuronList{k}.y;
                end
                    plot(coox,cooy,'.','Color',cmap(k,:),'Markersize',0.5*(1200/(5+abs( S - cooz))));
            end
            hold off;
        end
    end

% -=< Mouse scroll wheel callback function >=-
    function mouseScroll (object, eventdata,number)
        if nargin < 3
            UPDN = eventdata.VerticalScrollCount;
        else
            UPDN = number;
        end
        if (strcmp(View,'A'))
            S = S - UPDN;
        else
            S = S - 6*UPDN;
        end
        if (S < 1)
            S = 1;
        elseif (S > sno)
            S = sno;
        end
        if sno > 1
            set(shand,'Value',S);
            set(stxthand, 'String', sprintf('Slice# %d / %d',S, sno));
        else
            set(stxthand, 'String', '2D image');
        end
        
        if strcmp(View,'A')
            if viewport == 1
                slice = squeeze(Horizontal(:,:,S));
            else
                slice = squeeze(Sagittal(round(S*6),:,:))';
            end
        elseif strcmp(View,'S')
            if viewport == 1
                slice = squeeze(Horizontal(:,S,:))';
            else
                slice = squeeze(Sagittal(:,S,:));
            end
        else
            if viewport == 1
                slice = squeeze(Horizontal(S,:,:))';
            else
                slice = squeeze(Sagittal(:,:,round(S/6)));
            end
        end
        if size(slice,1) < 400
            temp = imresize(slice, [size(slice,1)*6.0 size(slice,2)]);
        elseif size(slice,2) < 400
            temp = imresize(slice, [size(slice,1) size(slice,2)*6.0]);
        else
            temp = slice;
        end
        
        K = get(gca,'children');
        set(K(end),'cdata',temp);
        delete(K(1:(end-1)));
        caxis([Rmin Rmax])
        if sno > 1
            set(stxthand, 'String', sprintf('Slice# %d / %d',S, sno));
        else
            set(stxthand, 'String', '2D image');
        end
        
        
        if blobs == 1
            cmap = lines(2560); hold on;
            for k = 1:numel(NeuronList)
                if strcmp(View,'A')
                    coox = NeuronList{k}.x; cooy = NeuronList{k}.y; cooz = NeuronList{k}.z; 
                elseif strcmp(View,'S')
                    coox = NeuronList{k}.y; cooy = NeuronList{k}.z*6; cooz = NeuronList{k}.x; 
                else
                    coox = NeuronList{k}.x; cooy = NeuronList{k}.z*6; cooz = NeuronList{k}.y;
                end
                    plot(coox,cooy,'.','Color',cmap(k,:),'Markersize',0.5*(1200/(5+abs( S - cooz))));
            end
            hold off;
        end
%         set(get(gca,'children'),'cdata',temp)
%         set(get(gca,'children'),'cdata',Horizontal(:,:,S))
    end

    function scroll_zoom(varargin)
        h = hittest(gcf);
        if isempty(h), return, end
        switch get(h,'Type')
          case 'axes'
            currAx = h;
          case 'image'
            currAx = get(h,'Parent');
          case 'line'
            currAx = get(h,'Parent'); 
          otherwise
            return
        end
%         if ~any(currAx == h_ax), return, end
        [x,y] = getAbsCoords(currAx);
        if ~coordsWithinLimits(currAx,x,y), return, end
        [x_rel, y_rel] = abs2relCoords(currAx, x, y);
        sc = varargin{2}.VerticalScrollCount;
        zoomFactor = abs(sc)*(1+20/100)^sign(sc);
        for i = currAx
          new_xlim_rel = ([0,1] - x_rel) * zoomFactor + x_rel;
          new_ylim_rel = ([0,1] - y_rel) * zoomFactor + y_rel;
          [new_xlim(1) new_ylim(1)] = rel2absCoords(currAx, new_xlim_rel(1), new_ylim_rel(1));
          [new_xlim(2) new_ylim(2)] = rel2absCoords(currAx, new_xlim_rel(2), new_ylim_rel(2));
          setNewLimits(currAx, new_xlim, new_ylim)
        end
    end


  function setNewLimits(ax, xlim, ylim)
    validX = ~any(isnan(xlim)) && ~any(isinf(xlim)) && diff(xlim)>0;
    if strcmp(get(ax,'XScale'),'log')
      validX = validX && ~isinf(xlim(2)/xlim(1));
    end
    if validX
      set(ax, 'Xlim', xlim);
    else
      if strcmp(tX.Running, 'off')
        old_color = get(ax, 'YColor');
        set(ax,'YColor','r');
        tX.TimerFcn = @(x,y)set(ax,'YColor',old_color);
        start(tX);
      end
    end
    
    validY = ~any(isnan(ylim)) && ~any(isinf(ylim)) && diff(ylim)>0;
    if strcmp(get(ax,'YScale'),'log')
      validY = validY && ~isinf(ylim(2)/ylim(1));
    end
    if validY
      set(ax, 'Ylim', ylim);
    else
      if strcmp(tY.Running, 'off')
        old_color = get(ax, 'XColor');
        set(ax,'XColor','r');
        tY.TimerFcn = @(x,y)set(ax,'XColor',old_color);
        start(tY);
      end
    end
  end

  function tf = coordsWithinLimits(h_ax, x, y)
    % check if the given point (x,y) is within the limits of the axis h_ax
    XLim = get(h_ax, 'xlim');
    YLim = get(h_ax, 'ylim');
    tf = x>XLim(1) && x<XLim(2) && y>YLim(1) && y<YLim(2);
  end

  function [x, y, z] = getAbsCoords(h_ax)
    crd = get(h_ax, 'CurrentPoint');
    x = crd(2,1);
    y = crd(2,2);
    z = crd(2,3);
  end

  function [x, y] = rel2absCoords(h_ax, x_rel, y_rel)
    XLim = get(h_ax, 'xlim');
    if strcmp(get(h_ax, 'XScale'), 'log')
      x = exp(x_rel*log(XLim(2)/XLim(1))+log(XLim(1)));
    else
      x = x_rel*diff(XLim)+XLim(1);
    end
    YLim = get(h_ax, 'ylim');
    if strcmp(get(h_ax, 'YScale'), 'log')
      y = exp(y_rel*log(YLim(2)/YLim(1))+log(YLim(1)));
    else
      y = y_rel*diff(YLim)+YLim(1);
    end
  end

  function [x_rel, y_rel] = abs2relCoords(h_ax, x, y)
    XLim = get(h_ax, 'xlim');
    if strcmp(get(h_ax, 'XScale'), 'log')
      x_rel = log(x/XLim(1))/log(XLim(2)/XLim(1));
    else
      x_rel = (x-XLim(1))/(XLim(2)-XLim(1));
    end
    YLim = get(h_ax, 'ylim');
    if strcmp(get(h_ax, 'YScale'), 'log')
      y_rel = log(y/YLim(1))/log(YLim(2)/YLim(1));
    else
      y_rel = (y-YLim(1))/(YLim(2)-YLim(1));
    end
  end

% -=< Mouse button released callback function >=-
    function mouseRelease (object,eventdata)
        MouseStat = get(gcbf, 'SelectionType');
        if strcmp(MouseStat,'alt')
            set(gcf, 'WindowButtonMotionFcn', '')
        end
    end

% -=< Mouse click callback function >=-
    function mouseClick (object, eventdata)
        MouseStat = get(gcbf, 'SelectionType');
        if (MouseStat(1) == 'a')        %   RIGHT CLICK
            InitialCoord = get(0,'PointerLocation');
            set(gcf, 'WindowButtonMotionFcn', @WinLevAdj);
        elseif (MouseStat(1) == 'e')        %   MIDDLE CLICK
            
        end
    end

% -=< Window and level mouse adjustment >=-
    function WinLevAdj(varargin)
        PosDiff = get(0,'PointerLocation') - InitialCoord;

        Win = Win + PosDiff(1) * WLAdjCoe * FineTuneC(get(ChBxhand,'Value')+1);
        LevV = LevV - PosDiff(2) * WLAdjCoe * FineTuneC(get(ChBxhand,'Value')+1);
        if (Win < 1)
            Win = 1;
        end

        [Rmin, Rmax] = WL2R(Win,LevV);
        caxis([Rmin, Rmax])
        set(lvalhand, 'String', sprintf('%6.0f',LevV));
        set(wvalhand, 'String', sprintf('%6.0f',Win));
        InitialCoord = get(0,'PointerLocation');
    end

% -=< Window and level text adjustment >=-
    function WinLevChanged(varargin)

        LevV = str2double(get(lvalhand, 'string'));
        Win = str2double(get(wvalhand, 'string'));
        if (Win < 1)
            Win = 1;
        end

        [Rmin, Rmax] = WL2R(Win,LevV);
        caxis([Rmin, Rmax])
    end

% -=< Window and level to range conversion >=-
    function [Rmn Rmx] = WL2R(W,L)
        Rmn = L - (W/2);
        Rmx = L + (W/2);
        if (Rmn >= Rmx)
            Rmx = Rmn + 1;
        end
    end

% -=< Window and level auto adjustment callback function >=-
    function AutoAdjust(object,eventdata)
        Win = double(max(Horizontal(:))-min(Horizontal(:)));
        Win (Win < 1) = 1;
        LevV = double(min(Horizontal(:)) + (Win/2));
        [Rmin, Rmax] = WL2R(Win,LevV);
        caxis([Rmin, Rmax])
        set(lvalhand, 'String', sprintf('%6.0f',LevV));
        set(wvalhand, 'String', sprintf('%6.0f',Win));
    end

% -=< load existing neuron list (strcat(outputfile,'.mat')) >=-
    function loadX(object,eventdata)
        temp = load(strcat(outputfile,'.mat'));
        NeuronList = temp.NeuronList;
    end

% -=< save neuron list ('NeuronList.mat'), overwrites without asking >=-
    function saveX(object,eventdata)
        save(strcat(outputfile,'.mat'),'NeuronList');
    end

% -=< Delete last neuron in list callback function >=-
    function deleteX(object,eventdata)
        if strcmp(questdlg('Do you really want to delete the last neuron in the list?'),'Yes');
            if numel(NeuronList) > 0;
                NeuronList(numel(NeuronList)) = [];
            end
        end
        mouseScroll(0,0,0);
    end

% -=< Create new neuron callback function >=-
    function createX(object,eventdata)
        [x_1,x_2] = ginput(1); 
        if View == 'C'
            x1 = x_1;
            x3 = x_2/6;
            x2 = S;
        elseif View == 'A'
            x1 = x_1;
            x2 = x_2;
            x3 = S;
        elseif View == 'S'
            x2 = x_1;
            x3 = x_2;
            x1 = S;
        end
        x = round(x1); y = round(x2); z = round(x3);
%         end
        index = numel(NeuronList)+1;
        NeuronList{index}.x = x;
        NeuronList{index}.y = y;
        NeuronList{index}.z = z;
        NeuronList{index}.x_abs = (x-772)*0.4233; % in micron from center
        NeuronList{index}.y_abs = (y-773)*0.4233; % in micron
        NeuronList{index}.z_abs = (z-55)*2.5; % in micron
        clockX = clock;
        NeuronList{index}.date = [date,',',32,num2str(clockX(4)),'h',num2str(clockX(5)),'min.'];
        NeuronList{index}.comment = get(Commhand,'string');
        
        disp(['Neuron number',32,num2str(index),',',32,'position',32,mat2str([NeuronList{index}.x NeuronList{index}.y NeuronList{index}.z]),32,'px, and',32,mat2str(round([NeuronList{index}.x_abs NeuronList{index}.y_abs NeuronList{index}.z_abs])),32,'micron from center']);
        
%         
%         [index NeuronList{index}.x NeuronList{index}.y NeuronList{index}.z]
    end

% -=< Create new neuron callback function (2) >=-
    function keypresser(object,eventdata)
        password = eventdata.Key;
        if strcmp(password,'x')
            createX();
            mouseScroll(0,0,0);
        elseif strcmp(password,'q')
            set (gcf, 'WindowScrollWheelFcn', @mouseScroll);
        elseif strcmp(password,'w')
            set (gcf, 'WindowScrollWheelFcn', @scroll_zoom);
        elseif strcmp(password,'b')
            blobs = abs(1-blobs);
            mouseScroll(0,0,0);
        elseif strcmp(password,'s')
            viewport = 1;
            mouseScroll(0,0,0);
        elseif strcmp(password,'d')
            viewport = 2;
            mouseScroll(0,0,0);
        end
    end

% -=< Axial view callback function >=-
    function AxialView(object,eventdata)
        if View == 'S'
            S_s = S;
        elseif View == 'C'
            S_c = S;
        end            
        View = 'A';
        
        S = S_a;
        sno = sno_a;
        cla(hdl_im);
        hdl_im = axes('position',[0,0.2,1,0.8]);
        
        if viewport == 1
            slice = squeeze(Horizontal(:,:,S));
        else
            slice = squeeze(Sagittal(S,:,:));
        end
        if size(slice,1) < 400
            temp = imresize(slice, [size(slice,1)*6.0 size(slice,2)]);
        elseif size(slice,2) < 400
            temp = imresize(slice, [size(slice,1) size(slice,2)*6.0]);
        else
            temp = slice;
        end        
        imshow(temp, [Rmin Rmax])

        if sno > 1
            shand = uicontrol('Style', 'slider','Min',1,'Max',sno,'Value',S,'SliderStep',[1/(sno-1) 10/(sno-1)],'Position', S_Pos,'Callback', {@SliceSlider, Horizontal});
            stxthand = uicontrol('Style', 'text','Position', Stxt_Pos,'String',sprintf('Slice# %d / %d',S, sno), 'BackgroundColor', [0.8 0.8 0.8], 'FontSize', SFntSz);
        else
            stxthand = uicontrol('Style', 'text','Position', Stxt_Pos,'String','2D image', 'BackgroundColor', [0.8 0.8 0.8], 'FontSize', SFntSz);
        end
        
        caxis([Rmin Rmax])
        if sno > 1
            set(stxthand, 'String', sprintf('Slice# %d / %d',S, sno));
        else
            set(stxthand, 'String', '2D image');
        end
        
        if size(slice,1) < 400
            temp = imresize(slice, [size(slice,1)*6.0 size(slice,2)]);
        elseif size(slice,2) < 400
            temp = imresize(slice, [size(slice,1) size(slice,2)*6.0]);
        else
            temp = slice;
        end  
            
        set(get(gca,'children'),'cdata',temp)
        set (gcf, 'ButtonDownFcn', @mouseClick);
        set(get(gca,'Children'),'ButtonDownFcn', @mouseClick);
        mouseScroll(0,0,0);
    end

% -=< Sagittal view callback function >=-
    function SagittalView(object,eventdata)
        if View == 'A'
            S_a = S;
        elseif View == 'C'
            S_c = S;
        end            
        View = 'S';
        S = S_s;
        sno = sno_s;
        if viewport == 1
            slice = squeeze(Horizontal(:,S,:));
        else
            slice = squeeze(Sagittal(:,S,:));
        end
        if size(slice,1) < 400
            temp = imresize(slice, [size(slice,1)*6.0 size(slice,2)]);
        elseif size(slice,2) < 400
            temp = imresize(slice, [size(slice,1) size(slice,2)*6.0]);
        else
            temp = slice;
        end
        
        cla(hdl_im);
        hdl_im = axes('position',[0,0.2,1,0.8]);
        
        imshow(temp, [Rmin Rmax])

        if sno > 1
            shand = uicontrol('Style', 'slider','Min',1,'Max',sno,'Value',S,'SliderStep',[1/(sno-1) 10/(sno-1)],'Position', S_Pos,'Callback', {@SliceSlider, Horizontal});
            stxthand = uicontrol('Style', 'text','Position', Stxt_Pos,'String',sprintf('Slice# %d / %d',S, sno), 'BackgroundColor', [0.8 0.8 0.8], 'FontSize', SFntSz);
        else
            stxthand = uicontrol('Style', 'text','Position', Stxt_Pos,'String','2D image', 'BackgroundColor', [0.8 0.8 0.8], 'FontSize', SFntSz);
        end
        
        caxis([Rmin Rmax])
        if sno > 1
            set(stxthand, 'String', sprintf('Slice# %d / %d',S, sno));
        else
            set(stxthand, 'String', '2D image');
        end

        if viewport == 1
            slice = squeeze(Horizontal(:,S,:));
        else
            slice = squeeze(Sagittal(:,S,:));
        end
        if size(slice,1) < 400
            temp = imresize(slice, [size(slice,1)*6.0 size(slice,2)]);
        elseif size(slice,2) < 400
            temp = imresize(slice, [size(slice,1) size(slice,2)*6.0]);
        else
            temp = slice;
        end
        set(get(gca,'children'),'cdata',temp)
        set (gcf, 'ButtonDownFcn', @mouseClick);
        set(get(gca,'Children'),'ButtonDownFcn', @mouseClick);
        mouseScroll(0,0,0);
    end

% -=< Coronal view callback function >=-
    function CoronalView(object,eventdata)
        if View == 'A'
            S_a = S;
        elseif View == 'S'
            S_s = S;
        end            
        View = 'C';
%         if viewport == 1
%             Horizontal = ImgCr;
%         else
%             Horizontal = Sagittal;
%         end
        S = S_c;
        sno = sno_c;
        cla(hdl_im);
        hdl_im = axes('position',[0,0.2,1,0.8]);
        if viewport == 1
            slice = squeeze(Horizontal(S,:,:))';
        else
            slice = squeeze(Sagittal(:,:,round(S/6)));
        end
        if size(slice,1) < 400
            temp = imresize(slice, [size(slice,1)*6.0 size(slice,2)]);
        elseif size(slice,2) < 400
            temp = imresize(slice, [size(slice,1) size(slice,2)*6.0]);
        else
            temp = slice;
        end
        

        imshow(temp, [Rmin Rmax])

        if sno > 1
            shand = uicontrol('Style', 'slider','Min',1,'Max',sno,'Value',S,'SliderStep',[1/(sno-1) 10/(sno-1)],'Position', S_Pos,'Callback', {@SliceSlider, Horizontal});
            stxthand = uicontrol('Style', 'text','Position', Stxt_Pos,'String',sprintf('Slice# %d / %d',S, sno), 'BackgroundColor', [0.8 0.8 0.8], 'FontSize', SFntSz);
        else
            stxthand = uicontrol('Style', 'text','Position', Stxt_Pos,'String','2D image', 'BackgroundColor', [0.8 0.8 0.8], 'FontSize', SFntSz);
        end
        
        caxis([Rmin Rmax])
        if sno > 1
            set(stxthand, 'String', sprintf('Slice# %d / %d',S, sno));
        else
            set(stxthand, 'String', '2D image');
        end
        
        if viewport == 1
            slice = squeeze(Horizontal(S,:,:))';
        else
            slice = squeeze(Sagittal(:,:,round(S/6)));
        end
        if size(slice,1) < 400
            temp = imresize(slice, [size(slice,1)*6.0 size(slice,2)]);
        elseif size(slice,2) < 400
            temp = imresize(slice, [size(slice,1) size(slice,2)*6.0]);
        else
            temp = slice;
        end
            
        set(get(gca,'children'),'cdata',temp)
        set (gcf, 'ButtonDownFcn', @mouseClick);
        set(get(gca,'Children'),'ButtonDownFcn', @mouseClick);
        mouseScroll(0,0,0);
    end

end