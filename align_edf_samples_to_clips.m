root_data_p = '/Volumes/external3/data/changlab/jamie/free-viewing';

data_p = fullfile( root_data_p, 'data' );
vid_p = fullfile( root_data_p, 'videos' );
preproc_p = fullfile( root_data_p, 'edf_samples' );
shared_utils.io.require_dir( preproc_p );

% sesh_dir = '08142023';
sesh_dirs = arrayfun( @(x) sprintf('08%d2023', x), 14:31, 'un', 0 );
% sesh_dirs = { '08142023' };

sesh_ps = fullfile( data_p, sesh_dirs );
sesh_ps = sesh_ps(cellfun(@(x) exist(x, 'file'), sesh_ps) > 0);
task_file_ps = find_task_data_files( sesh_ps );

parfor i = 1:numel(task_file_ps)
%%
  
fprintf( '\n %d of %d', i, numel(task_file_ps) );

task_file = shared_utils.io.fload( task_file_ps{i} );
sesh_p = fileparts( task_file_ps{i} );
edf_file = Edf2Mat( fullfile(sesh_p, task_file.edf_file_name) );

sync_info = extract_edf_sync_info( ...
  edf_file.Events.Messages, task_file.edf_sync_times );

dend_table = shared_utils.io.fload( fullfile(data_p, 'dendro_table.mat') );
clip_table = append_clip_meta_data_to_clip_table( ...
  task_file.params.target_clips, task_file.clip_table, dend_table );

clip_table.timestamp = NaT( rows(clip_table), 1 );
clip_table.timestamp(:) = task_file.time0_timestamp;

%%

edf_info = compute_edf_sample_traces( edf_file, sync_info, clip_table, vid_p );
clip_table.edf_info = edf_info;

save_p = fullfile( preproc_p, sprintf('edf_samples_%d.mat', i) );
do_save( save_p, clip_table );

end

%%

function do_save(preproc_p, clip_table)

save( preproc_p, 'clip_table' );

end

function edf_infos = compute_edf_sample_traces(edf_file, sync_info, clip_table, vid_p)

edf_infos = cell( height(clip_table), 1 );

for i = 1:height(clip_table)  
  %%
  
  fprintf( '\n\t %d of %d', i, height(clip_table) );
  
  vid_name = clip_table.video_filename{i};
  vid_reader = VideoReader( fullfile(vid_p, vid_name) );
  vid_fps = vid_reader.FrameRate;
  
  clip_index = sync_info{sync_info.video_time == clip_table.start(i), 'clip_index'};
  assert( numel(clip_index) == 1 && clip_index == i );
  match_clip = sync_info(sync_info.clip_index == clip_index, :);
  
  [~, start_ind] = min( match_clip.video_time );
  [~, stop_ind] = max( match_clip.video_time );  
  edf_start_t = match_clip.edf_time(start_ind);
  edf_stop_t = match_clip.edf_time(stop_ind);
  
  edf_t0_ind = find( edf_file.Samples.time == edf_start_t );
  edf_t1_ind = find( edf_file.Samples.time == edf_stop_t );
  
  %%
  
  vid_ts = nan( edf_t1_ind - edf_t0_ind + 1, 1 );
  vid_fis = nan( size(vid_ts) );
  
  edf_p = nan( size(vid_ts, 1), 2 );
  edf_ps = nan( size(edf_p, 2), 1 );
  
  for j = 1:edf_t1_ind-edf_t0_ind + 1
    ti = edf_t0_ind + j - 1;
    
    curr_edf_t = edf_file.Samples.time(ti);
    edf_px = edf_file.Samples.posX(ti);
    edf_py = edf_file.Samples.posY(ti);
    ps = edf_file.Samples.pupilSize(ti);
    
    vid_ts(j) = shared_utils.sync.cinterp( ...
      curr_edf_t, match_clip.edf_time, match_clip.video_time );
    vid_fis(j) = floor( vid_ts(j) * vid_fps );
    
    edf_p(j, :) = [ edf_px, edf_py ];
    edf_ps(j) = ps;
  end
  
  edf_infos{i} = struct( ...
      'position', edf_p, 'pupil_size', edf_ps ...
    , 'video_time', vid_ts ...
    , 'video_frame', vid_fis ...
  );
end

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