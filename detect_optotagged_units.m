function optotagged_units = detect_optotagged_units(filename,optotagged_threshold)
    % Lisez le fichier Excel
    data = readtable(filename);

    % Vérifiez les unités avec un reliability_scores >= optotagged_threshold
    selected_rows = data.reliability_scores >= optotagged_threshold;

    % Récupérez les numéros des unités optotaggées
    optotagged_units = data.units(selected_rows);
end
