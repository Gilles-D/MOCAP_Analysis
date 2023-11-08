%% Plots the Firing Locations of Neural Units
%  This function plots the firing locations of neural units based on their 
%  spike times and mouse position data. It uses Kernel Density Estimation 
%  (KDE) to plot the firing location for each unit. You can optionally pause 
%  between each plot.
%
% -------------------------------------------------------------------------
%% Syntax:
%   plot_firing_location(spike_times, time_axis, x_ref, y_ref, 
%                        observation_labels, spike_times_labels, 
%                        observations, p)
%
% -------------------------------------------------------------------------
%% Inputs:
%   spike_times(MATRIX):
%       A matrix where each column corresponds to the spike times of a
%       specific neural unit. NaN values should be filtered out.
%
%   time_axis(VECTOR):
%       A vector containing the time points at which observations were made.
%       Must be sorted in ascending order.
%
%   x_ref(STR):
%       The label used to compute a vector representing the x-coordinate
%       of the mouse at each time  point in 'time_axis'.
%
%   y_ref(STR):
%       The label used to compute a vector representing the y-coordinate
%       of the mouse at each time  point in 'time_axis'.
%
%   observation_labels(ANY TYPE) - Optional:
%       No implementation details are provided in the code. Placeholder for 
%       future features or debugging.
%
%   spike_times_labels(CELL ARRAY OF STRINGS) - Optional:
%       A cell array containing labels for each neural unit in 'spike_times'.
%
% 	observations(Matrix):
%   	Matrix containing behavior observations. Each row corresponds to an
%   	observation and each column to a behavior type
%
%   p(SCALAR) - Optional:
%       A scalar value specifying the pause duration in seconds between each 
%       unit's plot. Default is 0 (no pause).
%
% -------------------------------------------------------------------------
%% Outputs:
%   None
%       The function produces plots and does not return any value.
%
% -------------------------------------------------------------------------
%% Extra Notes:
%
% * The function internally handles NaN values in 'spike_times'.
% * The function uses 'interp1' to find the nearest indices in 'time_axis'
%   for each spike time.
% * 'kde' function is used for plotting the Kernel Density Estimation.
% * The function will pause for 'p' seconds after plotting each unit if 'p' 
%   is specified.
%
% -------------------------------------------------------------------------
%% Examples:
% * Plot with default pause time
%   plot_firing_location(spike_times, time_axis, x_ref, y_ref);
%
% * Plot with a 2-second pause between each unit
%   plot_firing_location(spike_times, time_axis, x_ref, y_ref, [], [], 2);
%
% -------------------------------------------------------------------------
%% Author(s):
%   Antoine Valera
%
% -------------------------------------------------------------------------
%                               Notice
%
% Notice Content will be added later. Leave a blank line here.
% -------------------------------------------------------------------------
% Revision Date:
%   20-09-2023
% -------------------------------------------------------------------------
% See also: 
%   kde, interp1
%
% TODO : 
% * Implement features or debugging options for 'observation_labels'.
% * Better handle 'spike_times_labels' for labeling plots.

function KDEs = plot_firing_location(spike_times, time_axis, x_ref, y_ref ,observation_labels, spike_times_labels, observations, p,save_path)  
    if nargin < 3 || isempty(x_ref)
        x_ref = 'back1_x';
    end
    if nargin < 4 || isempty(y_ref)
        y_ref = 'back1_y';
    end

    mouse_x = observations(:, find(strcmp(observation_labels, x_ref)));
    if isempty(mouse_x)
        mouse_x = observations(:, find(contains(observation_labels, x_ref),1));
    end
    mouse_y = observations(:, find(strcmp(observation_labels, y_ref)));
    if isempty(mouse_y)
        mouse_y = observations(:, find(contains(observation_labels, y_ref),1));
    end
    
    if nargin < 8 || isempty(p)
        p = 0                           ; 
    end
    
    %% Main loop to plot the firing location of each unit
    % Loop through all units in 'spike_times'. For each unit, compute the nearest
    % indices and plot KDE (Kernel Density Estimation).
    KDEs = {};
    for unit = 1:size(spike_times, 2)
        
        %% Data preparation for the current unit
        % - Filter out NaN values
        % - Remove time points that are outside of the time_axis range
        best_st = spike_times(:, unit)                   ;
        best_st = best_st(1:find(~isnan(best_st), 1, 'last')) ;
        best_st(best_st < time_axis(1)) = []             ;
        best_st(best_st > time_axis(end)) = []           ;
        
        %% Calculate nearest indices
        % Interpolate the time points to nearest indices in 'time_axis'
        indices         = 1:length(time_axis)                ;
        nearest_indices = interp1(time_axis, indices, best_st, 'nearest') ;
        nearest_indices(isnan(best_st)) = NaN               ;
        nearest_indices = nearest_indices(~isnan(nearest_indices)) ;
        
        %% KDE and plotting
        % Perform KDE and plot for the current unit
        KDEs{unit} = kde(mouse_x(nearest_indices), mouse_y(nearest_indices), true, [], [], ...
                    mouse_x, mouse_y, spike_times_labels{unit}, x_ref, y_ref) ;
                
        %% Save
            % Save figure
        if ~isempty(save_path)
            label = strrep(spike_times_labels{unit}, '\', '');
            filename = fullfile(save_path, ['KDE_', label, '.svg']);
            saveas(gcf, filename, 'svg');
        end
        
        %% Pause for visualization
        % Pause for 'p' seconds before plotting the next unit
        pause(p)                                        ;        
    end
end
