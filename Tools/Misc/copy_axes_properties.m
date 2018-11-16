function ax2 = copy_axes_properties(ax1,ax2)
%UNTITLED9 Summary of this function goes here
%   Detailed explanation goes here
ax2.XLim = ax1.XLim;
ax2.YLim = ax1.YLim;
ax2.ZLim = ax1.ZLim;
ax2.OuterPosition = ax1.OuterPosition;
ax2.Position = ax1.Position;
ax2.View = ax1.View;
ax2.PlotBoxAspectRatio = ax1.PlotBoxAspectRatio;
ax2.Color = 'none';
linkaxes([ax1,ax2]);
end

