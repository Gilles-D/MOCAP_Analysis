% function optotagged_units = detect_optotagged_units(filename,optotagged_threshold)
%     % Lisez le fichier Excel
%     data = readtable(filename);
% 
%     % Vérifiez les unités avec un reliability_scores >= optotagged_threshold
%     selected_rows = data.reliability_scores >= optotagged_threshold;
% 
%     % Récupérez les numéros des unités optotaggées
%     optotagged_units = data.units(selected_rows);
% end

function optotagged_indices = detect_optotagged_units(filename, zscore_threshold, success_rate_threshold)
    % Lisez le fichier Excel
    data = readtable(filename);

    % Vérifiez les unités avec un Zscore >= zscore_threshold et un pourcentage de succès >= success_rate_threshold
    selected_rows = (data.('Z_score') >= zscore_threshold) & (data.('x_Success') >= success_rate_threshold);


    % Récupérez les unités qui répondent aux critères
%     optotagged_units = data.units(selected_rows);
    optotagged_indices = find(selected_rows);
end

