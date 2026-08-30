function [marker_data, marker_map] = load_VICON_raw(filename)

    log_msg('load_VICON_raw', 'loading %s', filename);

    opts = detectImportOptions(filename);

    opts.VariableTypes(:)   = {'double'};
    opts.VariableNamesLine  = 3;
    opts.DataLines          = [6, Inf];
    opts.VariableNamingRule = 'preserve';

    T = readtable(filename, opts);
    
    T(:, 1:2) = [];

    vars = T.Properties.VariableNames;
    
    numMarkers  = floor(width(T) / 3);
    numFrames   = height(T);
    markerNames = cell(numMarkers, 1);
    marker_data = zeros(numMarkers, 3, numFrames);

    for i = 1:numMarkers
    
        col = (i - 1) * 3 + 1;
        markerNames{i} = matlab.lang.makeValidName(vars{col});
        marker_data(i, :, :) = T{:, col:col + 2}' / 1000;
    
    end

    % 19 mm markers
    marker_data(:, 3, :) = marker_data(:, 3, :) + 0.0095;

    marker_map = table((1:numMarkers)', markerNames, 'VariableNames', {'Index', 'Name'});

    log_msg('load_VICON_raw', '%d markers, %d raw frames', numMarkers, numFrames);
end
