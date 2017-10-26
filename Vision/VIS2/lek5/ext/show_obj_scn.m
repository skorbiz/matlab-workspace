function show_obj_scn(object, scene, varargin)
%SHOW_OBJ_SCN   Visualize a 3D object model and scene.
%
% SHOW_OBJ_SCN(object, scene) shows 3D object data by green color and 3D
% scene data by red color.
%
% SHOW_OBJ_SCN(object, scene, title) puts a title on top of the 3D plot.
%
% SHOW_OBJ_SCN(object, scene, title, h) puts a title on top of the 3D plot
% and shows the plot in an existing axis handle.
%
% Input:
%   object - The n-by-3 object point set
%   scene - The m-by-3 scene point set
%   title - The plot title
%   h - The axis handle to plot into
%
% Author: Anders Glent Buch

% Get optional input arguments which are [title axis_handle]
axis_title = '';
h = -1;
if nargin > 2
    axis_title = varargin{1};
    if nargin > 3
        h = varargin{2};
    end
end

% If no axis handle was specified, create a new figure
if h == -1
    figure;
    h = gca;
end

% Set h as current axis and add title
axes(h);
title(axis_title);

hold on

% Show object/scene
% scatter3(object(:,1), object(:,2), object(:,3), 10, 'go')
% scatter3(scene(:,1), scene(:,2), scene(:,3), 10, 'ro')
plot3(object(:,1), object(:,2), object(:,3), 'go')
plot3(scene(:,1), scene(:,2), scene(:,3), 'ro')
% Show origo
scatter3([0], [0], [0], 'k*');
line([0 0.1], [0 0], [0 0], 'Color', 'r', 'LineWidth', 2);
line([0 0], [0 0.1], [0 0], 'Color', 'g', 'LineWidth', 2);
line([0 0], [0 0], [0 0.1], 'Color', 'b', 'LineWidth', 2);
% Adjust axis and view
axis equal
view(3)
grid

hold off

end