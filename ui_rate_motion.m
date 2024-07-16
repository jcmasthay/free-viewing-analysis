function ratings = ui_rate_motion(movie_file)

if ( nargin < 1 )
  movie_file = '/Users/Nick/Downloads/jamie-efe-pose-track/GX010025.avi';
end

do_clean = onCleanup( @my_cleanup );

vr = VideoReader( movie_file );
num_frames = floor( vr.Duration * max(24, vr.FrameRate) );

f1 = figure(1);
clf( f1 );
ax1 = gca;

f2 = uifigure( 'name', 'slider' );
slider = uislider( 'Parent', f2 );

start_button = uibutton( 'Parent', f2, 'Position', [0, 0, 200, 50] );
start_button.ButtonPushedFcn = @enable_rating;
start_button.Text = 'begin';

quit_button = uibutton( 'Parent', f2, 'Position', [200, 0, 200, 50] );
quit_button.ButtonPushedFcn = @do_quit;
quit_button.Text = 'quit';

ratings = nan( num_frames, 1 );
enabled = false;
should_quit = false;

while ( ~enabled )
  pause( 1/60 );
end

fi = 1;
while ( ~ptb.util.is_esc_down && fi <= num_frames && ~should_quit )  
  frame = read( vr, fi );  
  ratings(fi) = slider.Value;
  
  cla( ax1 );
  imshow( frame, 'parent', ax1 );
  title( ax1, sprintf('Rating = %0.3f; Frame %d of %d' ...
    , ratings(fi), fi, num_frames) ); 
  drawnow;
  
  fi = fi + 1;
end

my_cleanup();

function do_quit(varargin)
  should_quit = true;
end

function enable_rating(varargin)
  enabled = true;
end

function my_cleanup()
  all_fig = findall(0, 'type', 'figure');
  delete( all_fig );
end

end