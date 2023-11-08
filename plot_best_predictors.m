%% Plot Best Predictors and Create Subplots
%   This function initializes two figures for plotting best and good
%   predictors of a given behavior, and populates them with subplots using
%   the `createSubplots` helper function. Each figure shows average rates
%   for good and best predictors for a selected behavior.
%
% -------------------------------------------------------------------------
%% Syntax:
%   plot_best_predictors(predictors, good_plus, good_minus, best_plus, best_minus, observations, subset, beh_idx, behaviour_subset)
%
% -------------------------------------------------------------------------
%% Inputs:
%   predictors(Matrix [Timepoints x Neurons]):
%       Matrix containing the predictor variables, usually the firing rates
%       of neurons.
%
%   good_plus(Vector [1 x n]):
%       Index vector specifying the columns in 'predictors' that are good
%       predictors with positive influence.
%
%   good_minus(Vector [1 x n]) - Optional:
%       Index vector specifying the columns in 'predictors' that are good
%       predictors with negative influence.
%
%   best_plus(Vector [1 x n]) - Optional:
%       Index vector specifying the columns in 'predictors' that are best
%       predictors with positive influence.
%
%   best_minus(Vector [1 x n]) - Optional:
%       Index vector specifying the columns in 'predictors' that are best
%       predictors with negative influence.
%
%   observations(Matrix [Timepoints x Features]):
%       Matrix containing the observed behavior variables.
%
%   subset(Vector [1 x m]):
%       Subset of behaviors to focus on, specified as indices.
%
%   beh_idx(Integer):
%       Index into 'behaviour_subset' specifying the behavior to plot.
%
%   behaviour_subset(Cell array of strings):
%       Cell array containing names of the behaviors that can be plotted.
%
% -------------------------------------------------------------------------
%% Outputs:
%   None. This function plots figures and does not return any outputs.
%
% -------------------------------------------------------------------------
%% Extra Notes:
%   * This function relies on the helper function `createSubplots` for
%     generating the subplots in each figure.
%   * Each subplot represents different aspects:
%       - First subplot shows the firing rate for positive predictors.
%       - Second subplot shows the specific behavior.
%       - Third subplot shows the firing rate for negative predictors.
%
% -------------------------------------------------------------------------
%% Examples:
% * Plot Best and Good Predictors
%   plot_best_predictors(pred_matrix, good_pos, good_neg, best_pos, best_neg, obs_matrix, [1, 2, 3], 1, {'Running', 'Eating'});
%
% -------------------------------------------------------------------------
%% Author(s):
%   Antoine Valera
%
% -------------------------------------------------------------------------
%                               Notice
%
%   Notice Content will be added later. Leave a blank line here.
% -------------------------------------------------------------------------
% Revision Date:
%   20-09-2023
%
% -------------------------------------------------------------------------
% See also: 

function plot_best_predictors(predictors, good_plus, good_minus, best_plus, best_minus, observations, subset, beh_idx, behaviour_subset)
    %% Initialize Figures
    % This block initializes two figures (663 and 664) and then calls the
    % createSubplots function to populate them with subplots.
    
    % Figure 663
    figure(663); clf(); hold on; set(gcf, 'color', 'w');
    createSubplots(predictors, good_plus,  good_minus, observations, subset, beh_idx, ...
                   ['Average rate for good predictors of ', behaviour_subset{beh_idx}]);

    % Figure 664
    figure(664); clf(); hold on; set(gcf, 'color', 'w');
    createSubplots(predictors, best_plus, best_minus, observations, subset, beh_idx, ...
                   ['Firing rate for best predictor of ',  behaviour_subset{beh_idx}]);
end

function createSubplots(predictors, good_plus, good_minus, observations, subset, beh_idx, tit)
    %% Create Subplots
    % This function creates three subplots and populates them with the relevant
    % data. The subplots are displayed in the current figure window.
    % - Subplot 311: Shows the firing rate for positive predictors
    % - Subplot 312: Shows the specific behavior
    % - Subplot 313: Shows the firing rate for negative predictors
    
    % Subplot 311
    subplot(311); 
    if all(good_plus <= size(predictors, 2))
        plot(mean(predictors(:, good_plus), 2), 'g');
    else
        error('Index in good_plus exceeds the number of columns in predictors.');
    end 
    set(gca, 'box', 'off'); 
    ylabel('F. rate  (+ \beta)');  % Explanatory comment
    
    % Subplot 312
    subplot(312); 
    plot(normalize(observations(:, subset(beh_idx))), 'k'); 
    set(gca, 'box', 'off'); 
    ylabel('behaviour');            % Explanatory comment
    
    % Subplot 313
    subplot(313); 
    plot(mean(predictors(:, good_minus), 2), 'r'); hold on; 
    set(gca, 'box', 'off'); 
    ylabel('F. rate  (- \beta)'); 
    xlabel('Timepoints');           % Explanatory comment
    
    % Set global title
    sgtitle(tit);
end
