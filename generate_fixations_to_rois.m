%%

vid_p = fullfile( fv_data_directory, 'videos' );
bbox_p = fullfile( fv_data_directory, 'detections' );
save_p = fullfile( fv_data_directory, 'fix_infos' );
shared_utils.io.require_dir( save_p );

samp_files = shared_utils.io.findmat( ...
  fullfile(fv_data_directory, 'edf_samples') ...
);

% samp_files = samp_files(1:20);

category_names = [ "animal", "human", "vehicle" ];

vid_infos = shared_utils.io.fload( fullfile(fv_data_directory, 'videos/vid_info.mat') );

%%

allow_overwrite = false;

%%

parfor si = 1:numel(samp_files)
  
%%

fprintf( '\n %d of %d', si, numel(samp_files) );

%%

fname = shared_utils.io.filenames( samp_files{si}, true );
dst_p = fullfile( save_p, fname );

if ( ~allow_overwrite && exist(dst_p, 'file') )
  continue
end
  
%%
  
clip_table = shared_utils.io.fload( samp_files{si} );
clip_table = convert_char_vars_to_string( clip_table );

% % @TODO: Haven't generated ROIs for these yet.
% has_mk = contains( clip_table.video_filename, 'Monkey Kingdom' );
% clip_table(has_mk, :) = [];

edf_infos = clip_table.edf_info;
vid_names = clip_table.video_filename;
% vid_names(:) = "Monkey Thieves S2E4.avi";

screen_dims = [1600, 900];
conf_threshold = 0.1;
category_types = { '1', '2', '3' };

%%

try  
  
detects_tbl = process_clip_table( ...
    edf_infos, vid_infos, vid_names, vid_p, bbox_p, screen_dims, conf_threshold ...
  , category_types ...
);

fprintf( '\n Saving %s ...', dst_p );
do_save( dst_p, detects_tbl );

catch err
  fprintf( '\n\n\n ||| %d failed', si );
  rethrow( err );
end

% per_file_outs{si} = tot_tbl;

end

%%

function do_save(p, tbl)
save( p, 'tbl' );
end

function detect_tbl = process_clip_table(...
  edf_infos, vid_infos, vid_names, vid_p, bbox_p, screen_dims, conf_threshold, category_types)

detect_tbl = table();

for i = 1:numel(edf_infos)
  %%
  
  fprintf( '\n\t %d of %d', i, numel(edf_infos) );
  
  bboxes = load( fullfile(bbox_p, sprintf('%s-bbox', vid_names(i)), 'all_bboxes.mat') );
%   vid_reader = VideoReader( fullfile(vid_p, vid_names(i)) );
  [~, loc] = ismember( vid_names(i), vid_infos.name );
  assert( loc > 0, 'No such video: "%s".', vid_names(i) );
  im_dims = vid_infos.size(loc, :);
  
%   im_dims = [ vid_reader.Width, vid_reader.Height ];
  
  %%
  
  fixs = edf_infos{i}.fixations;
  ib_detects = cell( numel(fixs), 1 );
  
  could_fix = zeros( numel(fixs), numel(category_types) );
  did_fix = zeros( size(could_fix) );
  
  dur_did_fix = zeros( numel(fixs), numel(category_types) );
  dur_could_fix = zeros( size(dur_did_fix) );
  
  weighted_dur_did_fix = zeros( size(dur_did_fix) );
  area_weight = zeros( size(dur_did_fix) );
  
  for j = 1:numel(fixs)
    fix = fixs(j);
    
    edf_px = fix.position(1);
    edf_py = fix.position(2);
    
    un_frames = unique( fix.video_frame(~isnan(fix.video_frame) & fix.video_frame > 0) );
    detects = bboxes.detections(un_frames);
    empties = cellfun(@(x) isequal(x, struct), detects);
    
    % remove empties
    detects(empties) = [];  
    un_frames(empties) = [];
    
    if ( ~isempty(detects) )
      un_frames = arrayfun( @(i) repmat(un_frames(i), 1, numel(detects{i})) ...
        , 1:numel(un_frames), 'un', 0 );
      detects = horzcat( detects{:} );
      un_frames = horzcat( un_frames{:} );
    end
    
    % remove low confidence
    confs = cellfun( @(x) x.conf, detects );
    detects(confs < conf_threshold) = [];
    un_frames(confs < conf_threshold) = [];
    
    detect_cats = cellfun( @(x) x.category, detects, 'un', 0 );
    
    rects = cellfun( ...
        @(x) bbox_to_pixel_rect(x.bbox, im_dims, screen_dims) ...
      , detects, 'un', 0 );
    
    is_ib = @(r) edf_px >= r(1) & edf_px < r(3) & edf_py >= r(2) & edf_py < r(4);
    ib_rect = cellfun( is_ib, rects );
    ib_detects(j) = { detects(ib_rect) };
    
    % for each unique frame, compute the area and area fraction for each 
    % category
    un_un = unique( un_frames );
    areas = zeros( size(un_un) );
    area_props = zeros( size(un_un) );
    
    for k = 1:numel(un_un)
      rs = rects(un_frames == un_un(k));
      tot_area = rect_area_calc( rs );
      areas(k) = tot_area;
      area_props(k) = tot_area ./ prod( screen_dims );
    end
    
    fix_areas = zeros( numel(category_types), numel(fix.video_frame) );
    fix_area_props = zeros( size(fix_areas) );
    for k = 1:numel(category_types)
      un_this_cat = unique( un_frames(strcmp(detect_cats, category_types{k})) );
      
      for h = 1:numel(un_this_cat)
        match_frame = fix.video_frame == un_this_cat(h);
        fix_areas(k, match_frame) = areas(h);
        fix_area_props(k, match_frame) = area_props(h);
      end
    end
    
    for k = 1:numel(category_types)      
      % indices of detections matching this category
      cat_inds = find( strcmp(detect_cats, category_types{k}) );
      
      tmp_could_fix = false( size(fix.video_frame) );
      tmp_did_fix = false( size(fix.video_frame) );
      
      for h = 1:numel(cat_inds)
        ci = cat_inds(h);
        match_frame = fix.video_frame == un_frames(ci);
        % there was a detection in this category to look to for these
        % (`match_frame)` frames.
        tmp_could_fix = tmp_could_fix | match_frame;
        
        % mask depending on whether the position was in bounds of this
        % detection.
        tmp_did_fix = tmp_did_fix | (match_frame * ib_rect(ci));
      end
      
      % proportionally smaller rois should have higher weight.
      weighted_dur_did_fix(j, k) = sum( double(tmp_did_fix(:)') ./ fix_area_props(k, :) );
      area_weight(j, k) = sum( 1 ./ fix_area_props(k, :) );
      
      dur_did_fix(j, k) = dur_did_fix(j, k) + sum( tmp_did_fix );
      dur_could_fix(j, k) = dur_could_fix(j, k) + sum( tmp_could_fix );
    end
    
    for k = 1:numel(category_types)
%       cellfun( @(x) x.category == category_types{k}, detects );
      matches_category = strcmp( detect_cats, category_types{k} );
      could_fix(j, k) = any( matches_category );
      did_fix(j, k) = any( ib_rect(matches_category) );
    end
  end
  
  src_index = repmat( i, numel(fixs), 1 );
  
  detect_tbl = [ detect_tbl; table( ...
    ib_detects, could_fix, did_fix ...
    , dur_could_fix, dur_did_fix ...
    , weighted_dur_did_fix ...
    , area_weight ...
    , src_index ...
  ) ];
  
end

end