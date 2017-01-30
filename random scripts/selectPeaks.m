function selectPeaks(~,event)
 
    global pos_manually
    clear keyword
    keyword = event.Key;
    switch(keyword)
        case 'x'
            [x,~] = ginput(1);
            if ~isempty(pos_manually)
                pos_manually(numel(pos_manually)+1) = x;
            else
                pos_manually(1) = x;
            end
        case 'q'
            pos_manually = pos_manually(1:end-1);
    end
end