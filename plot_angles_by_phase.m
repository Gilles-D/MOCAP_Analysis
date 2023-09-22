function plot_all_angles_by_phase(angles, phase)
    % Assumptions: 
    % angles is a vector of length N (here N = 15978)
    % phase is a vector of length N with values between -pi and pi

    % Initialize arrays to store angles for each phase
    phase_bins = -pi:0.1:pi;  % Modify the bin size if needed
    num_bins = length(phase_bins) - 1;
    angles_by_phase = cell(1, num_bins);
    
    % Group angles by phase
    for i = 1:num_bins
        idx = find(phase >= phase_bins(i) & phase < phase_bins(i+1));
        angles_by_phase{i} = angles(idx);
    end
    
    % Compute the mean angle for each phase bin
    mean_angles = zeros(1, num_bins);
    for i = 1:num_bins
        mean_angles(i) = mean(angles_by_phase{i});
    end
    
    % Plot all angles for each phase with transparency
    normalized_phase = linspace(0, 100, num_bins);  % Normalized phase in percentage
    hold on;
    for i = 1:num_bins
        plot(normalized_phase(i)*ones(size(angles_by_phase{i})), angles_by_phase{i}, '.r', 'MarkerSize', 8, 'Color', [1 0 0 0.1]);  % 0.1 is the alpha (transparency)
    end
    
    % Plot the mean angles on top
    plot(normalized_phase, mean_angles, 'o-', 'Color', [1 0 0]);
    xlabel('Normalized Phase (%)');
    ylabel('Angle');
    title('All Angles and Mean Angle by Phase');
    grid on;
    hold off;
end
