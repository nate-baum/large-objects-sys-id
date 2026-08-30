function p = get_marker(marker_data, marker_map, id)

    if isnumeric(id)
        idx = id;
    else
        idx = get_idx(marker_map, id);
    end

    p = squeeze(marker_data(idx, :, :));

end

function idx = get_idx(marker_map, marker_name)

    logical_idx = strcmp(marker_map.Name, marker_name);

    if ~any(logical_idx)
        error('Marker "%s" not found in the marker_map.', marker_name);
    end

    idx = marker_map.Index(logical_idx);

end
