function vid_info_table = generate_pose_tracking_video_info_table(csv_p)

if ( nargin < 1 )
  csv_p = '/Users/Nick/Downloads/jamie-efe-pose-track/Fellini Videos 2024.csv';
end

raw = readtable( csv_p );
dates = cellfun( @datestr, strrep(raw{1, :}, 'Pilot Day ', ''), 'un', 0 );
video_files = raw{2:end, :};

tot_files = {};
tot_dates = {};
tot_indices = [];
for i = 1:size(video_files, 2)
  first_empty = find( cellfun(@isempty, video_files(:, i)), 1 );
  if ( isempty(first_empty) )
    first_empty = size(video_files, 1); 
  else
    first_empty = first_empty - 1;
  end
  curr_files = video_files(1:first_empty, i);
  vid_index = 1:first_empty;
  
  tot_files = [ tot_files; curr_files ];
  tot_dates = [ tot_dates; repmat(dates(i), numel(curr_files), 1) ];
  tot_indices = [ tot_indices; vid_index(:) ];
end

tot_files = strrep( tot_files, 'O', '0' );
vid_ids = strrep( tot_files, '.MP4', '' );
vid_info_table = table( ...
  string(tot_files), string(tot_dates), string(vid_ids) ...
  , tot_indices, 'va', {'file', 'date', 'id', 'index'} );

end