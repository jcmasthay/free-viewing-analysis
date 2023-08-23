root_data_p = '/Volumes/external3/data/changlab/jamie/free-viewing';

data_p = fullfile( root_data_p, 'data' );
vid_p = fullfile( root_data_p, 'videos' );
bbox_p = fullfile( root_data_p, 'detections');

sesh_dir = '08142023';
task_files = find_task_data_files( fullfile(data_p, sesh_dir) );

task_file = shared_utils.io.fload( task_files{1} );
edf_file = shared_utils.io.fload( ...
  fullfile(data_p, sesh_dir, strrep(task_file.edf_file_name, '.edf', '.mat')) );

sync_info = extract_edf_sync_info( edf_file.Events.Messages, task_file.edf_sync_times );
ct = task_file.clip_table;

%%

accept_detect = @(x) x.conf >= 0.1;
screen_dims = [ task_file.window.Width, task_file.window.Height ];

store_detects = containers.Map();

look_durs = nan( height(ct), 1 );
for i = 1:height(ct)  
  fprintf( '\n %d of %d', i, height(ct) );
  
  vid_name = ct.video_filename{i};
  vid_reader = VideoReader( fullfile(vid_p, vid_name) );
  vid_fps = vid_reader.FrameRate;
  vid_dims = [ vid_reader.Width, vid_reader.Height ];
  
  clip_index = sync_info{sync_info.video_time == ct.start(i), 'clip_index'};
  assert( numel(clip_index) == 1 );
  match_clip = sync_info(sync_info.clip_index == clip_index, :);
  
%   match_clip = sync_info(sync_info.clip_index == ct.index(i), :);
  [vid_start_t, start_ind] = min( match_clip.video_time );
  [vid_stop_t, stop_ind] = max( match_clip.video_time );  
  edf_start_t = match_clip.edf_time(start_ind);
  edf_stop_t = match_clip.edf_time(stop_ind);
  
  edf_t0_ind = find( edf_file.Samples.time == edf_start_t );
  edf_t1_ind = find( edf_file.Samples.time == edf_stop_t );
  
  look_dur = 0;
  clip_t = vid_stop_t - vid_start_t;
  edf_t = (edf_t1_ind - edf_t0_ind + 1);
  
%   figure(1); clf;
%   subplot( 1, 2, 1 );
%   imshow( read(vid_reader, floor((vid_start_t + 3) * vid_fps)) );
%   subplot( 1, 2, 2 );
%   imshow( read(vid_reader, vid_stop_t * vid_fps) );
  
  for ti = edf_t0_ind:edf_t1_ind
    curr_edf_t = edf_file.Samples.time(ti);
    edf_px = edf_file.Samples.posX(ti);
    edf_py = edf_file.Samples.posY(ti);
    
    [~, nearest_ind] = min( abs(match_clip.edf_time - curr_edf_t) );
    vid_t = match_clip.video_time(nearest_ind);
    vid_fi = vid_t * vid_fps;
    
    bbox_filename = sprintf( 'bbox_%d.mat', vid_fi );
    bbox_file_p = fullfile( bbox_p, sprintf('%s-bbox', vid_name), bbox_filename );
    
    if ( isKey(store_detects, bbox_file_p) )
      bbox_file = store_detects(bbox_file_p);
      
    elseif ( ~exist(bbox_file_p, 'file') )
      continue
      
    else
      bbox_file = load( bbox_file_p );
      store_detects(bbox_file_p) = bbox_file;
    end
    
    detections = bbox_file.detections(cellfun(accept_detect, bbox_file.detections));
    detect_rects = cellfun( ...
      @(x) bbox_to_pixel_rect(x.bbox, vid_dims, screen_dims), detections, 'un', 0 );
    
    ib_rect = cellfun( @(r) edf_px >= r(1) & edf_px < r(3) & edf_py >= r(2) & edf_py < r(4) ...
      , detect_rects );
    look_dur = look_dur + any( ib_rect );
  end
  
  look_prop = look_dur / edf_t * clip_t;
  look_durs(i) = look_prop;
end

%%

figure(1);
clf;
hist( look_durs, 10 );
title( 'Duration of looking to animals' );

%%

function r = bbox_to_pixel_rect(bbox, im_dims, screen_dims)

p0 = bbox(1:2) .* im_dims;
wh = bbox(3:4) .* im_dims;

adj_x = (screen_dims(1) - im_dims(1)) * 0.5;
adj_y = (screen_dims(2) - im_dims(2)) * 0.5;

p0 = p0 + [ adj_x, adj_y ];
r = [ p0(1), p0(2), p0(1) + wh(1), p0(2) + wh(2) ];

end

function task_files = find_task_data_files(data_p)

% locate task data files in a directory and order them by earliest to
% latest in terms of data-collection time.

data_files = shared_utils.io.findmat( data_p );
file_names = shared_utils.io.filenames( data_files );

is_task_file = cellfun( @(x) numel(x) == 20, file_names );
task_dates = cellfun( @(x) strrep(x, '_', ':'), file_names, 'un', 0 );
task_nums = cellfun( @(x) datenum(datestr(x)), task_dates(is_task_file) );
[~, ord] = sort( task_nums );

task_files = data_files(is_task_file);
task_files = task_files(ord);

end