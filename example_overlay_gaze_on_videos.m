data_p = fullfile( fv_data_directory(), '08102023' );

mats = shared_utils.io.findmat( data_p );
f = load( mats{1} );
e = Edf2Mat( fullfile(data_p, f.file.edf_file_name) );
f = f.(char(fieldnames(f)));

clip_table = f.clip_table;
sync_info = extract_edf_sync_info( e.Events.Messages, f.edf_sync_times );

%%

begin = 1;
max_num_clips = 20;
vid_p = fullfile( project_directory, 'videos' );

win = ptb.Window( [0, 0, 1600, 900] );
win.SkipSyncTests = true;
open( win );

for i = begin:begin+min(max_num_clips, size(clip_table, 1))-1
  curr_clip = clip_table(i, :);
  vid_file_p = char( fullfile(vid_p, curr_clip.video_filename) );
  
  draw_cb = @(t) overlay_gaze( win, e, sync_info, i, t );
  play_movie( win, vid_file_p, curr_clip.start, curr_clip.stop, [], draw_cb );
end

close( win );

%%

function overlay_gaze(win, edf, sync_info, clip_index, vid_t)

match_clip = sync_info(sync_info.clip_index == clip_index, :);

if ( isempty(match_clip) || vid_t < min(match_clip.video_time) || vid_t > max(match_clip.video_time) )
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

x = (1 - f0) * x0 + f0 * x1;
y = (1 - f0) * y0 + f0 * y1;

s = 25;
r = [ x - s, y - s, x + s, y + s ];
Screen( 'FillOval', win.WindowHandle, [255, 0, 0], r );

end