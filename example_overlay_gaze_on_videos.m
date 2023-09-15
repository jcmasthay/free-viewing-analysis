block_index = 1;

data_p = 'C:\source\data\free_viewing\data\08102023';
bbox_p = 'C:\source\data\free_viewing\detections';

mats = shared_utils.io.findmat( data_p );
f = load( mats{block_index} );
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
conf_threshold = 0.6;

vid_p = fullfile( project_directory, 'videos' );
scram_vid_p = fullfile( vid_p, 'scrambled' );

win = ptb.Window( [0, 0, 1600, 900] );
win.SkipSyncTests = true;
open( win );

try

for i = begin:begin+min(max_num_clips, size(clip_table, 1))-1
  curr_clip = clip_table(i, :);

  if ( strcmp(curr_clip.block_type, 'C') )
    vid_file_p = char( fullfile(scram_vid_p, curr_clip.video_filename) );
  else
    vid_file_p = char( fullfile(vid_p, curr_clip.video_filename) );
  end

  vid_reader = VideoReader( vid_file_p );
  context.Value.video_reader = vid_reader;
  context.Value.bbox_p = fullfile( bbox_p, sprintf('%s-bbox', char(curr_clip.video_filename)) );
  context.Value.confidence_threshold = conf_threshold;

  draw_cb = @(t) overlay_gaze( win, context, e, sync_info, i, t, pupil_threshs );
  play_movie( win, vid_file_p, curr_clip.start, curr_clip.stop, [], draw_cb );
end

catch err
  warning( err.message );
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

vid_reader = context.Value.video_reader;
frame_index = floor( t0 * vid_reader.FrameRate );
im_dims = [ vid_reader.Width, vid_reader.Height ];

bbox_p = context.Value.bbox_p;
bbox_p = fullfile( bbox_p, sprintf('bbox_%d.mat', frame_index) );

if ( exist(bbox_p, 'file') )
  accept_detect = @(x) x.conf >= context.Value.confidence_threshold;
  bboxes = load( bbox_p );
  detections = bboxes.detections(cellfun(accept_detect, bboxes.detections));
  
  for i = 1:numel(detections)
    bbox = bbox_to_pixel_rect( detections{i}.bbox, im_dims, [win.Width, win.Height] );
    Screen( 'FrameRect', win.WindowHandle, [255, 255, 0], bbox, 4 );
  end
end

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
trail_w = 10;

for i = 1:size(pt, 1)-1
  frac_p = (i - 1) / max(1, size(pt, 1)-1);
  
  p0 = pt(i, 1:2);
  p1 = pt(i+1, 1:2);
  
  v = p1 - p0;
  vl = norm( v );
  v = v ./ vl;
  
  X = [-v(2); v(1)];
  Y = v(:);
  m = [ X, Y ];
  vs = (m * (verts .* [trail_w, vl])')' + p0 + Y' * vl * 0.5;  
  Screen( 'FillPoly', win.WindowHandle, [255 * frac_p, 0, 0], vs );
end

trail_len = 15;
update_trail( context.Value.position_trail, [x, y, psf], trail_len );

function r = make_rect(x, y, s, psf)
  adj_s = s * 0.5;
  s = s + adj_s * (psf * 2 - 1);
  r = make_sized_rect( x, y, s );
end

function r = make_sized_rect(x, y, s)
  r = [ x - s, y - s, x + s, y + s ];
end

end

function r = bbox_to_pixel_rect(bbox, im_dims, screen_dims)

p0 = bbox(1:2) .* im_dims;
wh = bbox(3:4) .* im_dims;

adj_x = (screen_dims(1) - im_dims(1)) * 0.5;
adj_y = (screen_dims(2) - im_dims(2)) * 0.5;

p0 = p0 + [ adj_x, adj_y ];
r = [ p0(1), p0(2), p0(1) + wh(1), p0(2) + wh(2) ];

end