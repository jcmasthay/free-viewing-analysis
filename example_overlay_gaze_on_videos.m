data_p = fullfile( fv_data_directory(), '08102023' );

mats = shared_utils.io.findmat( data_p );
f = load( mats{1} );
e = Edf2Mat( fullfile(data_p, f.file.edf_file_name) );
f = f.(char(fieldnames(f)));

clip_table = f.clip_table;
sync_info = extract_edf_sync_info( e.Events.Messages, f.edf_sync_times );
pupil_threshs = pupil_limits( e.Samples.pupilSize );

%%

position_trail = ptb.Reference( struct('history', []) );
context = ptb.Reference( struct('position_trail', position_trail) );

begin = 1;
max_num_clips = 20;
vid_p = fullfile( project_directory, 'videos' );

win = ptb.Window( [0, 0, 1600, 900] );
win.SkipSyncTests = true;
open( win );

for i = begin:begin+min(max_num_clips, size(clip_table, 1))-1
  curr_clip = clip_table(i, :);
  vid_file_p = char( fullfile(vid_p, curr_clip.video_filename) );
  
  draw_cb = @(t) overlay_gaze( win, context, e, sync_info, i, t, pupil_threshs );
  play_movie( win, vid_file_p, curr_clip.start, curr_clip.stop, [], draw_cb );
end

close( win );

%%

function threshs = pupil_limits(samples, props, num_bins)

if ( nargin < 3 )
  num_bins = 100;
end

if ( nargin < 2 )
  props = [ 0.125, 0.85 ];
end

[h, edges] = histcounts( samples, num_bins );
y = cumsum( h );
y = y ./ y(end);

ind0 = find( y > props(1), 1 );
ind1 = find( y > props(2), 1 );
threshs = [ edges(ind0), edges(ind1) ];

end

function update_trail(trail, curr_p, trail_len)

tv = trail.Value;

if ( size(tv.history, 1) < trail_len )
  tv.history(end+1, :) = curr_p;
else
  tv.history(1:end-1) = tv.history(2:end);
  tv.history(end, :) = curr_p;
end

trail.Value = tv;

end

function overlay_gaze(win, context, edf, sync_info, clip_index, vid_t, pupil_threshs)

match_clip = sync_info(sync_info.clip_index == clip_index, :);

if ( isempty(match_clip) || ...
     vid_t < min(match_clip.video_time) || ...
     vid_t > max(match_clip.video_time) )
  return
end

i0 = find( vid_t >= match_clip.video_time(1:end-1) & ...
           vid_t < match_clip.video_time(2:end), 1 );
i1 = min( i0 + 1, size(match_clip, 1) );

if ( isempty(i0) || isempty(i1) )
  return
end

t0 = match_clip.video_time(i0);
t1 = match_clip.video_time(i1);
f0 = (vid_t - t0) / (t1 - t0);

edf_t0_ind = edf.Samples.time == match_clip.edf_time(i0);
edf_t1_ind = edf.Samples.time == match_clip.edf_time(i1);

x0 = edf.Samples.posX(edf_t0_ind);
y0 = edf.Samples.posY(edf_t0_ind);
x1 = edf.Samples.posX(edf_t1_ind);
y1 = edf.Samples.posY(edf_t1_ind);

ps0 = edf.Samples.pupilSize(edf_t0_ind);
ps1 = edf.Samples.pupilSize(edf_t1_ind);

x = (1 - f0) * x0 + f0 * x1;
y = (1 - f0) * y0 + f0 * y1;
ps = (1 - f0) * ps0 + ps1 * y1;

ps = max( min(ps, pupil_threshs(2)), pupil_threshs(1) );
psf = (ps - pupil_threshs(1)) ./ (pupil_threshs(2) - pupil_threshs(1));
s = 25;

Screen( 'FillOval', win.WindowHandle, [255, 0, 0], make_rect(x, y, s, psf) );

% Draw history of gaze positions
pt = context.Value.position_trail.Value.history;
verts = [[-0.5, -0.5]; [0.5, -0.5]; [0.5, 0.5]; [-0.5, 0.5]];

for i = 1:size(pt, 1)-1
  p0 = pt(i, :);
  p1 = pt(i+1, :);
  
  v = p1(1:2) - p0(1:2);
  vl = norm( v );
  v = v ./ vl;
  theta = atan2( v(2), v(1) );
  m = rotation_2d( theta );
  vs = (m * (verts .* [ 5, vl ])')' + p0 + v * vl * 0.5;  
  frac_p = (i - 1) / max( 1, size(pt, 1) - 1 );
  
  Screen( 'FillPoly', win.WindowHandle, [255 * frac_p, 0, 0], vs );
  
%   Screen( 'FillOval', win.WindowHandle, [255 * frac_p, 0, 0] ...
%     , make_rect(prev_p(1), prev_p(2), s, prev_p(3)) );
end

trail_len = 10;
update_trail( context.Value.position_trail, [x, y, psf], trail_len );

function r = make_rect(x, y, s, psf)
  adj_s = s * 0.5;
  s = s + adj_s * (psf * 2 - 1);
  r = [ x - s, y - s, x + s, y + s ];    
end

end

function m = rotation_2d(theta, m3)

assert( numel(theta) == 1 );

m = [ cos(theta), -sin(theta)
      sin(theta), cos(theta) ];
    
if ( nargin > 1 && m3 )
  m3 = eye( 3 );
  m3(1:2, 1:2) = m;
  m = m3;
end

end