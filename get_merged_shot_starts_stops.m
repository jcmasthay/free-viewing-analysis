function [search_start_stop, merged] = get_merged_shot_starts_stops(samps)

% GET_MERGED_SHOT_STARTS_STOPS -- Get start and stop times of shots that
%   were or would be merged for survey data collection because the shot was 
%   less than 2 seconds long.
%
%   search_start_stop = get_merged_shot_starts_stops(samps) for the table
%   `samps` with variables 'start' and 'stop' returns an array 
%   `search_start_stop` of size `[rows(samps), 2]` whose first column 
%   contains the start times and second the stop times of the shots that 
%   would be shown to human subjects in the survey.

merge_dur = 2;
merge_each = {'block_type', 'Code', 'timestamp'};
merged = do_merge_shots( ...
    samps, 'start', 'stop' ...
  , findeach(samps, merge_each), merge_dur );

merged_search_start_stop = [ ceil(merged.start), floor(merged.stop) ];
search_start_stop = nan( rows(samps), 2 );

for i = 1:numel(merged.stop)
  assign_to = merged.merged_from{i};  
  search_start_stop(assign_to, :) = repmat( ...
    merged_search_start_stop(i, :), numel(assign_to), 1 );
end

if ( ~isempty(search_start_stop) && any(isnan(search_start_stop(end, :))) )
  pull_from = max( 1, rows(search_start_stop)-1 );
  search_start_stop(end, :) = search_start_stop(pull_from, :);
end

end