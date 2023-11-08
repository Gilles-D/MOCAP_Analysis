%% Run Machine Learning Algorithm on Mocap Data
%  This function is designed to run a machine learning algorithm on motion capture 
%  (mocap) data. It accepts standardized predictions and observations, along with
%  optional parameters for the time axis, observation labels, subsets, and machine 
%  learning parameters. The function performs iterative training and testing, then
%  visualizes the results with a bar chart.
%
% -------------------------------------------------------------------------
%% Syntax:
%  [results, mean_score, stats, ind_scores, subset, ml_params] = 
%       run_ml_mocap(   predictions, observations, time_axis,
%                       observation_labels, subset, ml_params)
%
% -------------------------------------------------------------------------
%% Inputs:
%  predictions(MATRIX):
%       Standardized predicted variables to be used in machine learning. It
%       is assumed to be standardized.
%
%  observations(MATRIX):
%       Standardized observed variables to be used in machine learning. It
%       is assumed to be standardized.
%
%  time_axis(MATRIX) - Optional:
%       Time axis data corresponding to the observations.
%
%  observation_labels(CELL ARRAY OF STRINGS) - Optional:
%       Labels corresponding to the observed variables. Used when filtering
%       by subset.
%
%  subset(EITHER MATRIX, STRING) - Optional:
%       * If subset is a string, then it filters the observations based on
%         labels in 'observation_labels'
%       * If subset is a matrix, then it should be a valid index array for
%         'observation_labels'
%
%  ml_params(STRUCT) - Optional:
%       Machine learning settings. Default settings are defined within the
%       function. For more details, see 'machine_learning_params'.
%
% -------------------------------------------------------------------------
%% Outputs:
%  results(CELL ARRAY):
%       Cell array storing results for each iteration of training and testing.
%
%  mean_score(DOUBLE):
%       Average score based on all iterations.
%
%  stats(STRUCT):
%       Statistics computed from the machine learning results.
%
%  ind_scores(MATRIX):
%       Individual scores for each iteration.
%
%  subset(ARRAY)
%       Array of behaviour indices that were selected for the analysis
%
%  ml_params(STRUCT):
%       Machine learning settings used for the analysis.
% -------------------------------------------------------------------------
%% Extra Notes:
%  * Observations and predictions should be standardized prior to being passed
%    into the function.
%  * The function executes the 'train_and_test' function iteratively based on
%    machine learning parameters defined in 'ml_params'.
%
% -------------------------------------------------------------------------
%% Examples:
% * Run with minimal inputs:
%   [results, mean_score] = run_ml_mocap(predictions, observations);
%
% * Run with custom subsets and ml_params:
%   [results, mean_score] = run_ml_mocap(predictions, observations, [], [], 'walking', custom_ml_params);
%
% -------------------------------------------------------------------------
%% Author(s):
%   Antoine Valera
%
% -------------------------------------------------------------------------
%                               Notice
%
%  Notice Content will be added later. Leave a blank line here
% -------------------------------------------------------------------------
% Revision Date:
%  19-09-2023
% -------------------------------------------------------------------------
% See also: 
%   train_and_test, bar_chart, machine_learning_params
%

function [results, mean_score, stats, ind_scores, subset, ml_params] = run_ml_mocap(predictions, observations, time_axis, observation_labels, subset, ml_params)
    % predictions and observations are expected to be standardized. Be careful when you pass the
    % input.

    %% Check for missing or empty subset argument and handle it accordingly
    % If 'subset' is not provided or is empty, default to using all columns of 'observations'.
    % If 'subset' is a string, then filter 'observations' based on the labels in 'observation_labels'.
    if nargin < 5 || isempty(subset)
        subset = 1:size(observations, 2)              ;
    elseif ischar(subset)
        subset = find(contains(observation_labels, subset));
    elseif any(subset > numel(observation_labels)) || any(subset <= 0)
        error('subsets must be a string (filtering observations based on the labels) or a valid indexing array of those labels');
    end

    %% Check for missing or empty ml_params argument and set default parameters
    % Prepare machine learning settings. For more details on options, see 'machine_learning_params'.
    if nargin < 6 || isempty(ml_params)
        ml_params = machine_learning_params('rendering'           , 1        , ...
                                            'holdout'             , 0.2      , ...
                                            'optimize_hyper'      , true     , ...
                                            'optimization_method' , 'manual' , ...
                                            'block_shuffling'     , 50       , ...
                                            'method'              , 'PLS' , ...
                                            'add_shuffle'         , false    , ...
                                            'is_shuffle'          , false    , ...
                                            'shuffling'           , 'behaviours' , ...
                                            'N_iter'              , 50);
    end

    %% Initialize variables
    predictor_idx = [] ;
    results       = {} ;
    %% Iterative training and testing
    % Execute 'train_and_test' function for each iteration and store results.
    for iter = 1:ml_params.N_iter
        fprintf(['Iteration : ', num2str(iter), '\n']);
        
        results{iter} = train_and_test(predictions'                 , ...
                                       observations(:, subset)'    , ...
                                       1:numel(time_axis)           , ...
                                       predictor_idx                , ...
                                       observation_labels(subset)   , ...
                                       observations(:, subset)'    , ...
                                       nanmedian(observations(:, subset), 2), ...
                                       ml_params);
    end

    %% Generate bar chart and collect statistics
    % Use 'bar_chart' function to visualize the results and collect statistical data.
    [mean_score, ~, fig_handle, stats, ind_scores] = bar_chart(results, '','','','', 1, true, ml_params);
end
