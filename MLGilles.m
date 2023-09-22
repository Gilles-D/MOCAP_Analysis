%% Load Excel file with behaviours and spike rates
data_filename = '\\equipe2-nas1\Public\DATA\Gilles\Spikesorting_August_2023\SI_Data\spikesorting_results\0022_01_08\kilosort3\curated\processing_data\0022_01_08_Mocap_Rates_obst_raw_norm_traj.xlsx';
[observations, predictors, time_axis, observation_labels, predictor_labels] = load_mocap_data(data_filename);

%% Load Excel file with spike times
spike_filename = '\\equipe2-nas1\Public\DATA\Gilles\Spikesorting_August_2023\SI_Data\spikesorting_results\0022_01_08\kilosort3\curated\processing_data\spike_times.xlsx';
[spike_times, spike_times_labels] = load_spike_times(spike_filename);

%% Quick Check for data integrity
if (numel(spike_times_labels) ~= numel(predictor_labels)) || ~all(ismember(spike_times_labels, predictor_labels))
    error('Unit list in spike times and spike rates do not match')
end

%% Get optotag infos
filename = '\\equipe2-nas1\Public\DATA\Gilles\Spikesorting_August_2023\SI_Data\spikesorting_results\0022_01_08\kilosort3\curated\processing_data/optotag_infos.xlsx';
optotagged_threshold = 60;
optotagged_units = detect_optotagged_units(filename,optotagged_threshold);
optotagged_labels = strcat('Unit\_', arrayfun(@num2str, optotagged_units, 'UniformOutput', false));
opto_indices = find(ismember(predictor_labels, optotagged_labels));


%% Interpolate data for small gaps based on gap duration
% small_gap_size_ms = 0.5;
% small_gap_size_tp = ceil(small_gap_size_ms / median(diff(time_axis)));
% for col = 1:size(observations,2)
%     observations(:, col) = fillmissing(observations(:, col),'pchip','EndValues', 'none','MaxGap',small_gap_size_tp);
% end

%% Remove spikes that outside the recorded behavioural data range
spike_times = trim_intertrial_events(time_axis, spike_times);


%% Fill observations gaps that span < 1/2 median step size
% QQ should be adjusted to compute a dynamic step size and fill gaps if
% less than 1/2 the current step size
[observations] = interpolate_gaps(observations, observation_labels, predictors, 'right_foot_x', 'back1_x');

%% Remove bad timepoints from observations and time axis
[observations, predictors, time_axis, spike_times] = ...
         trim_variables_for_missing_datapoints(observations, predictors, time_axis, spike_times);
     
%% Standardize data before machine learning
norm_obs = smoothdata(normalize(observations),'gaussian', 100);
norm_pred = smoothdata(normalize(predictors),'gaussian', 100);

% %% OPTIONAL Optotag subset
% norm_pred = norm_pred(:,opto_indices);

%% Start machine learning
[results, mean_score, stats, ind_scores, subset, ml_params] = run_ml_mocap(norm_pred, norm_obs, time_axis, observation_labels, '');


%% Get reference step cycle
reference_phase = get_step_cycle_phase_ref(observations, predictors, observation_labels, 'right_foot_x_norm');

%Get Stance state
ref_stance = 'right_foot_x'
swing_threshold = 50

stance_ref_idx = find(strcmp(observation_labels, ref_stance));
stance_ref_obs = observations(:,stance_ref_idx);
good_ref_obs = ~isnan(stance_ref_obs);
stance_ref_obs = stance_ref_obs(good_ref_obs);

stance_ref_speed = abs(diff(stance_ref_obs) / (1/200));
stance_indices = find(stance_ref_speed <= swing_threshold);
swing_indices = find(stance_ref_speed > swing_threshold);


% plot_angles_by_phase(observations(:,51),reference_phase);

% Figure 659: Plots the scores calculated from the model's coefficients, highlighting the mean score in black. Used for tuning assessment.
% Figure 663: Shows the mean of the predictors and observations for positively correlated neurons.
% Figure 664: Plots the predictors from the best positively and negatively correlated neurons.
% Figure 662: Shows scatter plots and linear fits for the mean firing rates of positively and negatively tuned neurons against the observed behavior.
% Figure 668: Histograms of mean angles and mean lengths.
% Figure 669: Shows the best-tuned cell based on the phase locking.
% Figure 777: KDE plots for each unit showing their firing behavior across positions.

%% Get phase of firing
behaviour_subset = observation_labels(subset);
sm_pred = smoothdata(predictors, 'gaussian', 100);

save_path = '\\equipe2-nas1\Public\DATA\Gilles\Spikesorting_August_2023\SI_Data\spikesorting_results\0022_01_08\kilosort3\curated\processing_data\plots\phaseresponse\'

single_unit_rendering = true;
[mean_angle, mean_magnitudes, n_events] = deal([]);        
for unit = 1:size(predictors, 2)
    [mean_angle(unit), mean_magnitudes(unit), ~, n_events(unit)] = plot_phase_of_reponse(observations(:, subset(1)), spike_times(:,unit), ['tuning for ',spike_times_labels{unit}], single_unit_rendering, reference_phase, time_axis,save_path,spike_times_labels{unit});
end

%%
save_path_scores = '\\equipe2-nas1\Public\DATA\Gilles\Spikesorting_August_2023\SI_Data\spikesorting_results\0022_01_08\kilosort3\curated\processing_data\plots\scores';

for beh_idx = 1:numel(behaviour_subset)
    if contains(behaviour_subset{beh_idx}, '')
        behaviour_subset{beh_idx} = strrep(behaviour_subset{beh_idx}, '_', '\_');
        score = [];
        for el = 1:ml_params.N_iter
            if strcmp(ml_params.method, 'linear')
                score(el, :) = results{el}.model{(beh_idx*2-1)}.Beta;
            elseif strcmp(ml_params.method, 'PLS')
                score(el, :) = results{el}.model{1}.Beta(:, beh_idx);
            end
        end
        
        [good_plus, best_plus, good_minus, best_minus, mean_score] = plot_score_figure(score, behaviour_subset, beh_idx, 50,save_path_scores,predictor_labels,optotagged_units);

        plot_best_predictors(predictors, good_plus, good_minus, best_plus, best_minus, observations, subset, beh_idx, behaviour_subset);
        plot_rate_vs_behaviour(observations, subset, beh_idx, behaviour_subset, predictors, good_plus, good_minus);

        min_n_event = 25;
        plot_cell_phase_locking(mean_angle, mean_magnitudes, n_events, min_n_event, predictors, observations, subset, beh_idx, behaviour_subset, predictor_labels);
    end
end


%% Get reference step cycle   
x_ref = 'right_foot_x';
y_ref = 'distance_from_obstacle_x';
save_path_KDE = '\\equipe2-nas1\Public\DATA\Gilles\Spikesorting_August_2023\SI_Data\spikesorting_results\0022_01_08\kilosort3\curated\processing_data\plots\KDE'
KDEs = plot_firing_location(spike_times, time_axis, x_ref, y_ref, observation_labels, spike_times_labels, observations, 1,save_path_KDE);
%%

mean_kde = cat(3, KDEs{:});
mean_kde = median(mean_kde, 3);

    mouse_x = observations(:, find(strcmp(observation_labels, x_ref)));
    if isempty(mouse_x)
        mouse_x = observations(:, find(contains(observation_labels, x_ref),1));
    end
    mouse_y = observations(:, find(strcmp(observation_labels, y_ref)));
    if isempty(mouse_y)
        mouse_y = observations(:, find(contains(observation_labels, y_ref),1));
    end
x_range = [min(mouse_x), max(mouse_x)] ;
y_range = [min(mouse_y), max(mouse_y)] ;

% Create a grid for the KDE
% Define min and max for x and y ranges and generate linearly spaced vectors
x_min = min(x_range)    ;
x_max = max(x_range)    ;
y_min = min(y_range)    ;
y_max = max(y_range)    ;


x_lin = linspace(x_min, x_max, 100) ;
y_lin = linspace(y_min, y_max, 100) ;
[X, Y] = meshgrid(x_lin, y_lin)     ;

% Flatten the grid matrices
X_flat = X(:) ;
Y_flat = Y(:) ;

figure()
contourf(X, Y, mean_kde, 'LineColor', 'none', 'Levels', 1e-5); hold on ;
% Plot all data points if provided
%scatter(mouse_x, mouse_y, 'Marker', 'o', 'MarkerFaceColor', 'w', 'MarkerFaceAlpha', 0.05, 'MarkerEdgeColor', 'none'); hold on ; % QQ: Consider setting up different styles for the scatter plot.
plot(mouse_x, mouse_y, '-','LineWidth',0.2,'Color',[0.5,0.5,0.5]); hold on ; % QQ: Consider setting up different styles for the scatter plot.

title('2D Kernel Density Estimation') ;
%set(gca, 'XDir', 'reverse');
xlabel(x_ref);
ylabel(y_ref);
colorbar ;
