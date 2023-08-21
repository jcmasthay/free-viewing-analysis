win = ptb.Window( [0, 0, 1280, 720] );
win.SkipSyncTests = true;
open( win );

start_times = [ 4 * 60 + 36 ];
stop_times = start_times + 15;

movie_file = fullfile( project_directory, 'videos/eg2.avi' );
play_movie(win, movie_file, start_times, stop_times, [], [])

close( win );

%%

transform_frame = @(im) imresize(im, [1280, 720]);
make_video_clips(movie_file, fullfile(fileparts(movie_file), 'upsized.avi') ...
  , [ start_times, stop_times ], 'Motion JPEG AVI', transform_frame );