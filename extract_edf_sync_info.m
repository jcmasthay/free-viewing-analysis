function info = extract_edf_sync_info(messages, mat_sync_times)

infos = messages.info(8:end);
info_ts = messages.time(8:end);

clip_indices = nan( numel(infos), 1 );
vid_timestamps = nan( size(clip_indices) );
sync_indices = nan( size(clip_indices) );

for i = 1:numel(infos)
  curr = strsplit( infos{i}, ' | ' );
  clip_indices(i) = str2double( curr{1} );
  vid_timestamps(i) = str2double( curr{2} );
  sync_indices(i) = str2double( curr{3} );
end

mat_time = mat_sync_times(sync_indices);

info = table( clip_indices, vid_timestamps, mat_time(:), sync_indices, info_ts(:) ...
  , 'va', {'clip_index', 'video_time', 'mat_time', 'mat_sync_index', 'edf_time'} );

end