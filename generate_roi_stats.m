%%

vid_p = fullfile( fv_data_directory, 'videos' );
bbox_p = fullfile( fv_data_directory, 'detections' );
samp_files = shared_utils.io.findmat( fullfile(fv_data_directory, 'edf_samples') );

%%

ct_vars = {'start', 'stop', 'video_filename', 'clip_id', 'affiliativeness', 'interactive_agency'};

tot_tbls = cell( numel(samp_files), 1 );

for i = 1:numel(samp_files)
  fprintf( '\n %d of %d', i, numel(samp_files) );
  samp = shared_utils.io.fload( samp_files{i} );
  samp_frames = cellfun( @(x) unique(x.video_frame), samp.edf_info, 'un', 0 );
  tot_tbls{i} = [ samp(:, ct_vars), table(samp_frames, 'va', {'sample_frames'}) ];
end

tot_tbl = vertcat( tot_tbls{:} );

[un_clips, ~, ic] = unique( tot_tbl(:, ct_vars) );
ic = groupi( ic );

samp_frames = cell( size(ic) );
for i = 1:numel(ic)
  samp_frames{i} = unique( cate1(tot_tbl.sample_frames(ic{i})) );
end
un_clips.sample_frames = samp_frames;

%%

conf_threshold = 0.1;
screen_dims = [1600, 900];

poss_cats = { '1', '2', '3' };
category_names = [ "animal", "human", "vehicle" ];

vid_names = string( un_clips.video_filename );
[vid_I, vid_C] = findeachv( vid_names );

roi_stat_tbl = table();

for i = 1:numel(vid_C)
  %%
  
  fprintf( '\n %d of %d', i, numel(vid_C) );  
  bboxes = load( fullfile(bbox_p, sprintf('%s-bbox', vid_C(i)), 'all_bboxes.mat') );
  vid_reader = VideoReader( fullfile(vid_p, vid_C(i)) );
  im_dims = [ vid_reader.Width, vid_reader.Height ];
  
  %%
  
  vi = vid_I{i};
  
  areas = zeros( numel(vi), numel(poss_cats) );
  mean_areas = zeros( size(areas) );
  freqs = zeros( size(areas) );
  
  for j = 1:numel(vi)
    frames = un_clips.sample_frames{vi(j)};
    detects = bboxes.detections(frames);
    
    detect_frames = arrayfun( @(i) repmat(frames(i), 1, numel(detects{i})) ...
        , 1:numel(frames), 'un', 0 );
    
    miss = cellfun( @(x) isequal(x, struct), detects );
    
    detects(miss) = [];
    detect_frames(miss) = [];
    
    detects = [ detects{:} ];
    detect_frames = [ detect_frames{:} ];
    
    confs = cellfun( @(x) x.conf, detects );
    detects(confs < conf_threshold) = [];    
    detect_frames(confs < conf_threshold) = [];
    rects = eachcell( @(x) bbox_to_pixel_rect(x.bbox, im_dims, screen_dims), detects );
    
    % compute areas and frequencies per frame and category
    [g_I, g_C] = findeach( [cellfun(@(x) str2double(x.category), detects(:)), detect_frames(:)], 1:2 );
    tmp_areas = cellfun( @(x) rect_area_calc(rects(x)) ./ prod(screen_dims), g_I );
    tmp_freqs = cellfun( @numel, g_I );
    
    % collapse within category
    [g_I, g_C] = findeach( g_C, 1 );
    areas_per_frame = cellfun( @(x) sum(tmp_areas(x)) / numel(frames), g_I );
    freqs_per_frame = cellfun( @(x) sum(tmp_freqs(x)) / numel(frames), g_I );
    
    [~, locs] = ismember( string(g_C), poss_cats );
    mean_areas(j, locs) = areas_per_frame;
    freqs(j, locs) = freqs_per_frame;
  end
  
  mean_areas(~isfinite(mean_areas)) = 0;
  areas(~isfinite(areas)) = 0;
  
  base_tbl = un_clips(vi, :);
  base_tbl.areas = areas;
  base_tbl.mean_areas = mean_areas;
  base_tbl.roi_freqs = freqs;
  base_tbl.frame_rate = repmat( vid_reader.FrameRate, size(vi) );
  
  roi_stat_tbl = [ roi_stat_tbl; base_tbl ];
end

%%

ind_sets = repmat( 1:rows(roi_stat_tbl), 1, numel(category_names) );
ind_ind = repelem( 1:numel(category_names), rows(roi_stat_tbl) );

per_cat_tbl = roi_stat_tbl(ind_sets, :);
per_cat_tbl.category = columnize( category_names(ind_ind) );

flat_vars = { 'areas', 'mean_areas', 'roi_freqs' };

for j = 1:numel(flat_vars)
  per_cat_tbl.(flat_vars{j}) = columnize( roi_stat_tbl{:, flat_vars(j)} );
end

%%

is_area = true;

if ( is_area )
  plt_vec = per_cat_tbl.mean_areas;
else
  plt_vec = per_cat_tbl.roi_freqs;

  if ( 0 )
    plt_vec = plt_vec ./ ((per_cat_tbl.stop - per_cat_tbl.start) .* per_cat_tbl.frame_rate);
  end
end

mask = rowmask( plt_vec );

[I, id, C] = rowsets( 4, per_cat_tbl ...
  , {'affiliativeness'}, {'interactive_agency'}, {'category'}, {} ...
  , 'mask', mask, 'to_string', true );
C = strrep( C, '_', ' ');

figure(2); clf;
[axs, hs, xs] = plots.simplest_barsets( plt_vec, I, id, C ...
  , 'summary_func', @nanmean ...
  , 'error_func', @plotlabeled.nansem ...
  , 'add_points', false ...
);

set( axs, 'xticklabelrotation', 10 );
shared_utils.plot.match_ylims( axs );

if ( is_area )
  ylabel( axs(1), 'Mean area proportion per frame' );
%   ylim( axs, [0, 1] );
else
  % ylabel( axs(1), 'Mean duration with roi type' );
  % ylabel( axs(1), 'Mean number of rois of type, per clip' );
  ylabel( axs(1), 'Mean number of rois per frame' );
end