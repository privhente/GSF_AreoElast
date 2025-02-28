function spy3D(S,param2,param3)
%SPY3D Visualize sparsity pattern with 3D bars.
%   SPY3D(S) plots the sparsity pattern of the matrix S in a
%   three-dimensional (3D) domain with matrix entries represented as bars
%   that are colored according to their value (height) and using a
%   colorbar.
%
%   See also SPY, BAR3 COLORBAR.
%
%   For more information, see <a href="matlab:
%   web('https://github.com/fjramireg/MatGen')">the MatGen Web site</a>.

%   Written by Francisco Javier Ramirez-Gil, fjramireg@gmail.com
%   Universidad Nacional de Colombia - Medellin
%   Created: 30/11/2018. Modified: 21/01/2019. Version: 1.3
%   updates:
%   --------
%   C. Hente, 21/02/2023
facealpha = 1.0;
dcolorbar = 5;

if exist('facealpha','var'); facealpha = param2; end
if exist('dcolorbar','var'); dcolorbar = param3; end

b = bar3(S);
cbh = colorbar;
for k = 1:length(b)
    zdata = b(k).ZData;
    for j = 1:size(zdata,1)/6
        max_vals = max(max(zdata(6*(j-1)+1:6*(j-1)+6,:)));
        min_vals = min(min(zdata(6*(j-1)+1:6*(j-1)+6,:)));
        abs_max_vals = max([abs(max_vals),abs(min_vals)]);
        si = 1;
        if abs(min_vals) >= abs(max_vals); si = -1; end
        zdata(6*(j-1)+1:6*(j-1)+6,:) = abs_max_vals*si;
    end
    b(k).CData = zdata;
    b(k).FaceColor = 'flat';
    b(k).FaceAlpha = facealpha;
end
cbh.Ticks = linspace(cbh.Limits(1),cbh.Limits(2),dcolorbar);
return