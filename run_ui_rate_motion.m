% determine which video to load
movie_filepath = '/path/to/movie.mp4';
% movie_filepath = '/Users/Nick/Downloads/jamie-efe-pose-track/GX010025.avi';

% determine where and whether to save
save_path = '~/Downloads';
do_save = true;

% run the ratings script
ratings = ui_rate_motion( movie_filepath );

% put the results into a table
[~, movie_id] = fileparts( movie_filepath );
ratings = table( ratings, strings(numel(ratings), 1) + movie_id ...
  , 'VariableNames', {'rating', 'movie_id'} );

% save the output
if ( do_save )
  save( fullfile(save_path, sprintf('ratings_%s.mat', movie_id)), 'ratings' );
end