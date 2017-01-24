function add_horizontal_axes( TTick, TTickLim, newAxisTime, newAxisValue, newAxisLabel, axisNo)
%add_horizontal_axes.m Add additional XTickLabels below the existing time
%axis
% Input
% TTick       : Time Tick values which ought to be labelled in MATLAB units
% TTickLim    : The min and max value of time axis in MATLAB units
% newAxisTime : [1xN] Time array of the new axis tick values
% newAxisValue: [1xN] Array of new axis values (Example: s/c position)
% newAxisLabel: The label of the new axis
% axisNo      : Integer values specifying the number of the axes. This
%               integer will determine the distance of the labels from 
%               the axis line
    
    DateNumBeg = TTickLim(1);
    DateNumEnd = TTickLim(2);
    normalizedTTick=(TTick-DateNumBeg)/(DateNumEnd-DateNumBeg);
    
    for iTick=1:1:length(TTick)
        iState(iTick)=find_time(newAxisTime, datestr(TTick(iTick)));
    end;
    
    delta = -0.11*(axisNo)-0.03*(axisNo-1); % 0.11 is the first axis tick label's vertical distance
                                             %0.03 is the vertical distance between
                                             %two horizontal axis tick labels
    
    for iTick=1:1:length(TTick)
        thisTick = normalizedTTick(iTick);
        txtHandle = text(thisTick, delta, num2str(round(newAxisValue(iState(iTick)),1)),'Units','normalized');
        set(txtHandle, 'HorizontalAlignment','center','VerticalAlignment','middle');
    end;
    
    txtHandle = text(-0.05, delta, newAxisLabel,'Units','normalized','FontSize',7);
    set(txtHandle, 'HorizontalAlignment','right','VerticalAlignment','middle');
    % -0.05 specifies the horizontal starting point of the right-aligned
    % axis label
end

