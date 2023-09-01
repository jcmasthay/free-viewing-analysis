root_data_p = fv_data_directory();

data_p = fullfile( root_data_p, 'data' );
vid_p = fullfile( root_data_p, 'videos' );
bbox_p = fullfile( root_data_p, 'detections');
preproc_p = fullfile( root_data_p, 'looking_proportions' );
shared_utils.io.require_dir( preproc_p );

% sesh_dir = '08142023';
sesh_dirs = arrayfun( @(x) sprintf('08%d2023', x), 14:31, 'un', 0 );
% sesh_dirs = { '08142023' };

sesh_ps = fullfile( data_p, sesh_dirs );
sesh_ps = sesh_ps(cellfun(@(x) exist(x, 'file'), sesh_ps) > 0);
task_file_ps = find_task_data_files( sesh_ps );

% task_file_ps = task_file_ps(1:8);

allow_overwrite = false;

parfor i = 1:numel(task_file_ps)
%%
  
fprintf( '\n %d of %d', i, numel(task_file_ps) );

task_file = shared_utils.io.fload( task_file_ps{i} );
sesh_p = fileparts( task_file_ps{i} );
edf_file = Edf2Mat( fullfile(sesh_p, task_file.edf_file_name) );
fname = shared_utils.io.filenames( task_file_ps{i} );

save_p = fullfile( preproc_p, sprintf('%s.mat', fname) );
if ( exist(save_p, 'file') && ~allow_overwrite )
  continue
end

sync_info = extract_edf_sync_info( ...
  edf_file.Events.Messages, task_file.edf_sync_times );

dend_table = shared_utils.io.fload( fullfile(data_p, 'dendro_table.mat') );
clip_table = append_clip_meta_data_to_clip_table( ...
  task_file.params.target_clips, task_file.clip_table, dend_table );

clip_table.timestamp = NaT( rows(clip_table), 1 );
clip_table.timestamp(:) = task_file.time0_timestamp;

%%

success = true;
try
  curr_look_prop_ct = compute_roi_looking_proportions( ...
    task_file, edf_file, sync_info, clip_table, vid_p, bbox_p );
catch err
  warning( err.message );
  success = false;
end

if ( ~success )
  continue
end

do_save( save_p, curr_look_prop_ct );

end

%%

function do_save(preproc_p, look_prop_ct)

save( preproc_p, 'look_prop_ct' );

end

function clip_table = compute_roi_looking_proportions(...
  task_file, edf_file, sync_info, clip_table, vid_p, bbox_p)

accept_detect = @(x) x.conf >= 0.1;
screen_dims = [ task_file.window.Width, task_file.window.Height ];

store_detects = containers.Map();

look_props = nan( height(clip_table), 1 );

for i = 1:height(clip_table)  
  %%
  
  fprintf( '\n\t %d of %d', i, height(clip_table) );
  
  vid_name = clip_table.video_filename{i};
  vid_reader = VideoReader( fullfile(vid_p, vid_name) );
  vid_fps = vid_reader.FrameRate;
  vid_dims = [ vid_reader.Width, vid_reader.Height ];
  
  clip_index = sync_info{sync_info.video_time == clip_table.start(i), 'clip_index'};
  assert( numel(clip_index) == 1 && clip_index == i );
  match_clip = sync_info(sync_info.clip_index == clip_index, :);
  
  [~, start_ind] = min( match_clip.video_time );
  [~, stop_ind] = max( match_clip.video_time );  
  edf_start_t = match_clip.edf_time(start_ind);
  edf_stop_t = match_clip.edf_time(stop_ind);
  
  edf_t0_ind = find( edf_file.Samples.time == edf_start_t );
  edf_t1_ind = find( edf_file.Samples.time == edf_stop_t );
  
  look_dur = 0;
  roi_dur = 0;
  
  %%
  
  for ti = edf_t0_ind:edf_t1_ind
    curr_edf_t = edf_file.Samples.time(ti);
    edf_px = edf_file.Samples.posX(ti);
    edf_py = edf_file.Samples.posY(ti);
    
    vid_t = shared_utils.sync.cinterp( ...
      curr_edf_t, match_clip.edf_time, match_clip.video_time );
    vid_fi = floor( vid_t * vid_fps );
    
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
    
    is_ib = @(r) edf_px >= r(1) & edf_px < r(3) & edf_py >= r(2) & edf_py < r(4);
    did_look = any( cellfun(is_ib, detect_rects) );
    
    look_dur = look_dur + did_look;
    % weight `did_look` by whether there was anything to look to.
    roi_dur = roi_dur + double( ~isempty(detections) );
  end
  
  %%
  
  if ( roi_dur ~= 0 )
    look_props(i) = look_dur / roi_dur;
  end
end

clip_table.look_props = look_props(:);

end

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