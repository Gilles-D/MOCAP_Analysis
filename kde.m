%% Performs 2D Kernel Density Estimation (KDE)
% 	This function performs a 2D Kernel Density Estimation (KDE) on the input data 
%   sets x_data and y_data. The function allows optional arguments for rendering 
%   the KDE plot and defining x and y ranges.
%
% -------------------------------------------------------------------------
%% Syntax:
% 	out = kde(x_data, y_data, rendering, x_range, y_range, all_x_data, all_y_data)
%
% -------------------------------------------------------------------------
%% Inputs:
% 	x_data(Vector) Double:
%   	Vector containing the x-coordinate data points.
%
% 	y_data(Vector) Double:
%   	Vector containing the y-coordinate data points.
%
% 	rendering(Scalar) Logical - Optional:
%       Flag to enable or disable rendering of the KDE plot.
%           * If true, the KDE plot will be rendered.
%           * If false, the KDE plot will not be rendered (default).
%
% 	x_range(Vector) Double - Optional:
%       Vector specifying the minimum and maximum x-values for the KDE.
%           * Default is based on the range of x_data.
%
% 	y_range(Vector) Double - Optional:
%       Vector specifying the minimum and maximum y-values for the KDE.
%           * Default is based on the range of y_data.
%
%   all_x_data(Vector) Double - Optional:
%       Vector containing all x-coordinate data points to be plotted as scatter.
%
%   all_y_data(Vector) Double - Optional:
%       Vector containing all y-coordinate data points to be plotted as scatter.
%
%   x_ref(STR):
%       The label used to compute a vector representing the x-coordinate
%       of the mouse at each time  point in 'time_axis'.
%
%   y_ref(STR):
%       The label used to compute a vector representing the y-coordinate
%       of the mouse at each time  point in 'time_axis'.
% -------------------------------------------------------------------------
%% Outputs:
% 	out(Matrix) Double:
%       Matrix representing the KDE density on the 2D grid defined by x_range 
%       and y_range.
%
% -------------------------------------------------------------------------
%% Extra Notes:
%
% * This function utilizes the 'ksdensity' function from MATLAB's Statistics
%   and Machine Learning Toolbox.
% * If the rendering flag is set to true, a contour plot is generated. 
%   Optionally, all data points can also be plotted if all_x_data and all_y_data 
%   are provided.
%
% -------------------------------------------------------------------------
%% Examples:
% * Basic usage
% 	out = kde(x_data, y_data);
%
% * With rendering
% 	out = kde(x_data, y_data, true);
%
% * With custom x and y range
% 	out = kde(x_data, y_data, false, [-5, 5], [-5, 5]);
%
% -------------------------------------------------------------------------
%% Author(s):
%   Antoine Valera
%
% -------------------------------------------------------------------------
%                               Notice
%
% -------------------------------------------------------------------------
% Revision Date:
% 	19-09-2023
% -------------------------------------------------------------------------
% See also: 
%	ksdensity

function out = kde(x_data, y_data, rendering, x_range, y_range, all_x_data, all_y_data, unit_label, x_ref, y_ref)
    %% Validate input arguments
    % Check if rendering flag, x_range and y_range are provided. If not, set defaults.
    if nargin < 3 || isempty(rendering)
        rendering = false                        ;
    end
    if nargin < 4 || isempty(x_range)
        x_range = [min(all_x_data), max(all_x_data)] ;
    end
    if nargin < 5 || isempty(y_range)
        y_range = [min(all_y_data), max(all_y_data)] ;
    end

    %% Create a grid for the KDE
    % Define min and max for x and y ranges and generate linearly spaced vectors
    x_min = min(x_range)    ;
    x_max = max(x_range)    ;
    y_min = min(y_range)    ;
    y_max = max(y_range)    ;

    x_lin = linspace(x_min, x_max, 100) ;
    y_lin = linspace(y_min, y_max, 100) ;
    [X, Y] = meshgrid(x_lin, y_lin)     ;

    %% Flatten the grid matrices
    X_flat = X(:) ;
    Y_flat = Y(:) ;

    %% Perform the KDE on the data
    % Utilize the 'ksdensity' function to compute the density
    [pdf_val, xi] = ksdensity([x_data(:), y_data(:)], [X_flat, Y_flat]);%, 'Bandwidth',5) ;

    %% Reshape the density to match the grid size
    Z = reshape(pdf_val, length(y_lin), length(x_lin)) ;

    %% Optionally render the KDE
    % Create a plot if rendering is enabled
    if rendering
        figure(777); clf() ; set(gcf, 'Color', 'w');
        
        subplot(5,1,[1:4])
        plot_fig(X,Y,Z,x_data,y_data,all_x_data,all_y_data, x_ref, y_ref, true);
        
        
        subplot(5,1,5)
        plot_fig(X,Y,Z,x_data,y_data,all_x_data,all_y_data, x_ref, y_ref, false); axis equal
        sgtitle(['Firing location for ',unit_label]);
    end

    %% Output the KDE density matrix
    out = Z ;
end

function plot_fig(X,Y,Z,x_data,y_data,all_x_data,all_y_data, x_ref, y_ref, plot_events)
    contourf(X, Y, Z, 'LineColor', 'none', 'Levels', 1e-5); hold on ;
    % Plot all data points if provided
    if nargin > 5
        scatter(all_x_data, all_y_data, 'Marker', 'o', 'MarkerFaceColor', 'w', 'MarkerFaceAlpha', 0.05, 'MarkerEdgeColor', 'none'); hold on ; % QQ: Consider setting up different styles for the scatter plot.
    end

    if plot_events
        scatter(x_data, y_data, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'none'); hold on ;
    end
    title('2D Kernel Density Estimation') ;
    %set(gca, 'XDir', 'reverse');
    xlabel(x_ref);
    ylabel(y_ref);
    colorbar ;
end