%% Plot Model Beta Weights and Identify Good and Best Predictors
% 	This function plots the Beta weights of a model and identifies both good 
%   and best predictors based on a given threshold. It also calculates the mean 
%   score from the input score data.
%
% -------------------------------------------------------------------------
%% Syntax:
% 	[good_plus, best_plus, good_minus, best_minus, mean_score] = 
%      plot_score_figure(score, behaviour_subset, beh_idx, good_predictors_cutoff)
%
% -------------------------------------------------------------------------
%% Inputs:
% 	score(DOUBLE):
% 		Array containing the score (Beta weights) for each predictor.
%
% 	behaviour_subset(CELL):
% 		Cell array containing behavior subset names.
%
% 	beh_idx(INTEGER):
% 		Index for the desired behavior subset within the 'behaviour_subset' array.
%
% 	good_predictors_cutoff(DOUBLE) - Optional:
% 		Threshold for identifying good predictors, expressed in percentage of 
%       max positive or min negative Beta value. Default is 50.
%
% -------------------------------------------------------------------------
%% Outputs:
% 	good_plus(VECTOR) DOUBLE:
% 		Indices of predictors with positive Beta values that exceed the threshold.
%
% 	best_plus(SCALAR) DOUBLE:
% 		Index of the best predictor with a positive Beta value.
%
% 	good_minus(VECTOR) DOUBLE:
% 		Indices of predictors with negative Beta values that exceed the threshold.
%
% 	best_minus(SCALAR) DOUBLE:
% 		Index of the best predictor with a negative Beta value.
%
%   mean_score(VECTOR) DOUBLE:
%       Mean score calculated from the input score array.
%
% -------------------------------------------------------------------------
%% Extra Notes:
% 	* If good_predictors_cutoff is provided as a value less than 1, a warning
%     will be shown to ensure that the user did not intend to enter it as a 
%     fraction.
%   * The function also plots the Beta weights, with individual predictors in 
%     gray and the mean in black.
%
% -------------------------------------------------------------------------
%% Examples:
% * Basic usage
% 	[good_plus, best_plus, good_minus, best_minus, mean_score] = 
%     plot_score_figure(score, behaviour_subset, 1);
%
% * With custom good_predictors_cutoff
% 	[good_plus, best_plus, good_minus, best_minus, mean_score] = 
%     plot_score_figure(score, behaviour_subset, 1, 60);
%
% -------------------------------------------------------------------------
%% Author(s):
% 	Antoine Valera
%
% -------------------------------------------------------------------------
%                               Notice
%  % paste license here
% -------------------------------------------------------------------------
% Revision Date:
% 	20-09-2023
% -------------------------------------------------------------------------
% See also: 


function [good_plus, best_plus, good_minus, best_minus, mean_score] = plot_score_figure(score, behaviour_subset, beh_idx, good_predictors_cutoff, save_path, predictor_labels,optotagged_units)
    %% Check input arguments and set defaults
    if nargin < 4 || isempty(good_predictors_cutoff)
        good_predictors_cutoff = 50;  % Predictors with Beta values above good_predictors_cutoff ...
                                      % ... % of the max are conserved (and split between + and - weights)
    elseif good_predictors_cutoff < 0 || good_predictors_cutoff > 100
        error('good_predictors_cutoff is a percentage of max response, and must be set between 0 and 100'); 
    elseif good_predictors_cutoff < 1
        warning('good_predictors_cutoff is expressed in %, not in fraction. You put a value < 1, make sure this is what you want');
    end

    %% Calculate mean score and standard error of the mean (SEM) from input score data
    mean_score = nanmean(score);
    sem = nanstd(score) / sqrt(size(score, 1)); % Standard Error of the Mean

    %% Calculate good and best predictors based on mean score
    good_plus  = find(mean_score > (max(mean_score(mean_score > 0)) * (good_predictors_cutoff / 100)));
    [~, best_plus ] = max(mean_score);  % Best predictor with positive Beta value

    good_minus = find(mean_score < (min(mean_score(mean_score < 0)) * (good_predictors_cutoff / 100)));
    [~, best_minus] = min(mean_score);  % Best predictor with negative Beta value

    %% Plot individual scores and mean score
    figure(659); clf();
    plot(score', 'Color', [0.8, 0.8, 0.8]); hold on;
    plot(mean_score, 'k', 'LineWidth', 2);
    ylim([-1, 1]);  % Verrouiller l'échelle des ordonnées entre -1 et 1
    set(gcf, 'color', 'w'); % Set figure background to white
    set(gca, 'box', 'off'); % Remove axis box
    xlabel('Unit'); ylabel('Model Beta Weight');
    title(strrep(strrep(behaviour_subset{beh_idx}, '\', ''),'_', ' ')); % Use behavior subset as the plot title
    xticks(1:length(predictor_labels));
    xticklabels(predictor_labels);
    xtickangle(45); % Rotate the labels for better visibility

    % Save the figure if save_path is provided
    if nargin >= 5 && ~isempty(save_path)
        label = strrep(behaviour_subset{beh_idx}, '\', '');
        saveas(gcf, fullfile(save_path, ['Score_' label '.svg']));
    end

    %% Plot the mean scores as bar chart with error bars
    optotagged_labels = strcat('Unit\_', arrayfun(@num2str, optotagged_units, 'UniformOutput', false));
    opto_indices = find(ismember(predictor_labels, optotagged_labels));

    figure(6591); clf();
    bar(mean_score, 'FaceColor', [0.8, 0.8, 0.8]); hold on;
    errorbar(1:length(mean_score), mean_score, sem, 'k.', 'LineWidth', 1.5);
    ylim([-1, 1]);  % Verrouiller l'échelle des ordonnées entre -1 et 1

    % Add asterisks for optotagged units
    marker_height = mean_score(opto_indices) + sem(opto_indices) + max(sem) * 0.1; % A bit above the error bar
    text(opto_indices, marker_height, '*', 'Color', 'b', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 12);

    set(gcf, 'color', 'w'); % Set figure background to white
    set(gca, 'box', 'off'); % Remove axis box
    xlabel('Unit'); ylabel('Model Beta Weight');
    title(strrep(strrep(behaviour_subset{beh_idx}, '\', ''),'_', ' ')); % Use behavior subset as the plot title
    xticks(1:length(predictor_labels));
    xticklabels(predictor_labels);
    xtickangle(45); % Rotate the labels for better visibility

    % Save the figure if save_path is provided
    if nargin >= 5 && ~isempty(save_path)
        label = strrep(behaviour_subset{beh_idx}, '\', '');
        saveas(gcf, fullfile(save_path, ['ScoreHisto_' label '.svg']));
    end

end


