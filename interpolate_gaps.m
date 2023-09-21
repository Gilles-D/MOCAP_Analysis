%% Interpolates gaps in mouse movement data based on reference variables
%  This function aims to interpolate gaps in the observational data of mouse
%  movements. It interpolates gaps in the data that are smaller than half
%  of a median step distance. It uses a step reference variable and a
%  progress reference variable to compute the median step size. This helps
%  improve data consistency when working with movement tracking data.
%
% -------------------------------------------------------------------------
%% Syntax:
% 	[observations] = interpolate_gaps(observations, observation_labels, 
%                                      predictors, ref_step_variable, 
%                                      ref_progress_variable)
%
% -------------------------------------------------------------------------
%% Inputs:
% 	observations(NxM Double Matrix):
%       Matrix containing observational data where N is the number of samples
%       and M is the number of features (observational variables).
%
% 	observation_labels(1xM Cell Array of Strings):
%       Array of strings that labels each feature in the observations matrix.
%
% 	predictors(NxP Double Matrix):
%       Matrix of predictor variables where N is the number of samples and
%       P is the number of predictor features.
%
% 	ref_step_variable(Optional String) - Optional:
%       Name of the variable that will serve as the reference step variable.
%       Default is 'right_foot_x'.
%
% 	ref_progress_variable(Optional String) - Optional:
%       Name of the variable that will serve as the reference for tracking
%       progress along the X-axis. Default is 'back1_x'.
%
% -------------------------------------------------------------------------
%% Outputs:
% 	observations(NxM Double Matrix):
%       Matrix containing interpolated observational data.
%
% -------------------------------------------------------------------------
%% Extra Notes:
% * It's important that the observation and predictor matrices are aligned
%   in terms of their row count (i.e., number of samples).
% * Gaps that are too large (more than half the median step size) are not
%   interpolated and a warning message is printed.
%
% -------------------------------------------------------------------------
%% Examples:
% * Interpolate gaps with default reference variables
% 	observations = interpolate_gaps(observations, observation_labels, 
%                                   predictors);
%
% * Interpolate gaps with custom reference variables
% 	observations = interpolate_gaps(observations, observation_labels, 
%                                   predictors, 'custom_step_ref', 
%                                   'custom_progress_ref');
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
% 	19-09-2023
% -------------------------------------------------------------------------
% See also: 
% 	plot_phase_of_reponse, filloutliers, fillmissing, findpeaks
%

function [observations] = interpolate_gaps(observations, observation_labels, predictors, ref_step_variable, ref_progress_variable)
    if nargin < 4 || isempty(ref_step_variable)
        ref_step_variable = 'right_foot_x';
    end
    if nargin < 5 || isempty(ref_progress_variable)
        ref_progress_variable = 'back1_x';
    end
    
    %% Interpolate gaps in data that last less than 1/2 a median step distance
    % Prolonged stationary periods end up losing the marker tracking, but as the mouse
    % didn't move, interpolation is relatively safe.
    
    %% Compute median step size
    % First, we find the median step size based on the specified reference step
    % variable (default is right_foot_x) and reference progress variable 
    % (default is back1_x).
    
    %% Validate array orientations
    % Make sure that observations are column vectors.
    if size(observations, 1) < size(observations, 2)
        observations = observations'                                           ; % Transpose if needed
    end

    %% Find reference step variable column
    %step_ref = find(contains(observation_labels, ref_step_variable))           ; % Find column matching reference step variable
    step_ref = find(strcmp(observation_labels, ref_step_variable))              ; % Find column matching exactly reference step variable
    step_ref = observations(:, step_ref)                                       ; % Extract that column
    good_step_ref = ~isnan(step_ref)                                           ; % Find non-NaN entries
    step_ref = step_ref(good_step_ref)                                         ; % Keep only non-NaN entries
    [~, ~, phase] = plot_phase_of_reponse(step_ref, predictors(:, 1), '', false); % Calculate phase of response
    [~, pk_loc] = findpeaks(phase, 'MinPeakProminence', pi)                    ; % Find peaks in phase to identify steps
    
    %% Find reference progress variable column
    x_progress_ref = find(contains(observation_labels, ref_progress_variable)) ; % Find column matching reference progress variable
    x_progress_ref = observations(:, x_progress_ref) * -1                      ; % Extract and negate that column
    x_progress_ref = diff(x_progress_ref)                                      ; % Compute the difference to get progress steps
    
    %% Fill outliers and missing data for x_progress_ref
    x_progress_ref = filloutliers(x_progress_ref, NaN)                          ; % Fill outliers with NaN
    x_progress_ref = fillmissing(x_progress_ref, 'pchip')                       ; % Interpolate missing data
    
    %% Calculate median step size
    x_progress_ref = cumsum(x_progress_ref)                                     ; % Cumulative sum to get total progress
    median_step = median(diff(x_progress_ref(pk_loc)))                           ; % Calculate median step size based on progress peaks
    x_progress_ref = [0; x_progress_ref]                                        ; % Prepend zero to align with original data size
    
    %% Try to interpolate data for small gaps
    for col = 1:size(observations, 2)
        
        %% Prepare current observation column
        current_obs = observations(:, col)                                      ; % Extract current column for observation
        current_obs = [0; current_obs; 0]                                       ; % Prepend and append zero to facilitate edge-case handling
        
        %% Identify NaN blocks
        isnanArr = isnan(current_obs)                                           ; % Identify NaN entries
        diffIsnan = diff(isnanArr)                                              ; % Compute difference to find starts and ends of NaN blocks
        
        %% Locate the start and end indices of each NaN block
        startPoints = find(diffIsnan == 1)                                      ; % Find starting points of NaN blocks
        endPoints = find(diffIsnan == -1) - 1                                   ; % Find ending points of NaN blocks
        
        %% Evaluate each NaN block for interpolation
        for i = 1:length(startPoints)
            
            %% Compute x-distance for the current NaN block
            delta_x_for_bout = x_progress_ref(endPoints(i)) - x_progress_ref(startPoints(i));
            
            %% Interpolate if x-distance is less than half of the median step size
            if delta_x_for_bout < median_step / 2
                end_win = min(endPoints(i) + 1, size(observations, 1))          ; % End window for interpolation
                observations(1:end_win, col) = fillmissing(observations(1:end_win, col), 'pchip', 'EndValues', 'none');
            else
                fprintf('NaN bout from index %d and ends at index %d is too large to be interpolated\n', startPoints(i), endPoints(i));
            end
        end
    end
end
