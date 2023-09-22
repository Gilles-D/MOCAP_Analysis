function plot_angles_by_phase(angles, phase)
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
    
    % Plot the mean angles against the normalized phase
    normalized_phase = linspace(0, 100, num_bins);  % Normalized phase in percentage
    plot(normalized_phase, mean_angles);
    xlabel('Normalized Phase (%)');
    ylabel('Mean Angle');
    title('Mean Angle by Phase');
    grid on;
end
