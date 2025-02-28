% +------------------------------------------------------+
% |      Grid Lines for 3-D, 4-D, 5-D and 6-D Plots      | 
% |              with MATLAB Implementation              |               
% |                                                      |
% | Author: Ph.D. Eng. Hristo Zhivomirov        11/25/19 | 
% +------------------------------------------------------+
%
% function: grid3(xgrid, ygrid, zgrid)
%
% Input:
% xgrid - x-grid tick locations on the x-axis;
% xgrid - y-grid tick locations on the y-axis;
% xgrid - z-grid tick locations on the z-axis.
% 
% Output:
% N/A 
%
% Note: the inputs may be vectors or 3D matrices 
% obatained, for instance, by the 'meshgrid' function.

function grid3(xgrid, ygrid, zgrid)

% input validation
xgrid = unique(xgrid(:));
ygrid = unique(ygrid(:));
zgrid = unique(zgrid(:)); 
validateattributes(xgrid, {'single', 'double'}, ...
                          {'real', 'nonnan', 'nonempty', 'finite'}, ...
                          '', 'xgrid', 1)
validateattributes(ygrid, {'single', 'double'}, ...
                          {'real', 'nonnan', 'nonempty', 'finite'}, ...
                          '', 'ygrid', 2)
validateattributes(zgrid, {'single', 'double'}, ...
                          {'real', 'nonnan', 'nonempty', 'finite'}, ...
                          '', 'zgrid', 3)

% remove all grid lines from the current axes
% grid off

% set the axis limits to the range of the data
% axis tight

% set the hold state to 'on'
% hold on

% determine the grids' limits
xgridmin = min(xgrid); xdatamax = max(xgrid);
ygridmin = min(ygrid); ydatamax = max(ygrid);
zgridmin = min(zgrid); zdatamax = max(zgrid);

% form the x-grid and the y-grid
for z = 1:length(zgrid)
    % x-grid
    for x = 1:length(xgrid)
        plot3([xgrid(x) xgrid(x)], ...
              [ygridmin ydatamax], ...
              [zgrid(z) zgrid(z)], '--', 'Color', [.5 .5 .5],'linewidth',0.1);
    end
    
    % y-grid
    for y = 1:length(ygrid)
        plot3([xgridmin xdatamax], ...
              [ygrid(y) ygrid(y)], ...
              [zgrid(z) zgrid(z)], '--', 'Color', [.5 .5 .5],'linewidth',0.1);
    end
end

% form the z-grid
for y = 1:length(ygrid)
    for x = 1:length(xgrid)
        plot3([xgrid(x) xgrid(x)], ...
              [ygrid(y) ygrid(y)], ...
              [zgridmin zdatamax], '--', 'Color', [.5 .5 .5],'linewidth',0.1);
    end
end

% set the hold state to 'off'
% hold off