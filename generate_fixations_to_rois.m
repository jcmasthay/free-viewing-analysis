%%

vid_p = fullfile( fv_data_directory, 'videos' );
bbox_p = fullfile( fv_data_directory, 'detections' );

samp_files = shared_utils.io.findmat( ...
  fullfile(fv_data_directory, 'edf_samples') ...
);

category_names = [ "animal", "human", "vehicle" ];

%%

per_file_outs = cell( numel(samp_files), 1 );

parfor si = 1:numel(samp_files)
  
%%

fprintf( '\n %d of %d', si, numel(samp_files) );
  
%%
  
clip_table = shared_utils.io.fload( samp_files{si} );
clip_table = convert_char_vars_to_string( clip_table );

% @TODO: Haven't generated ROIs for these yet.
has_mk = contains( clip_table.video_filename, 'Monkey Kingdom' );
has_e6 = contains( clip_table.video_filename, 'S2E6' );
clip_table(has_mk | has_e6, :) = [];

edf_infos = clip_table.edf_info;
vid_names = clip_table.video_filename;
% vid_names(:) = "Monkey Thieves S2E4.avi";

screen_dims = [1600, 900];
conf_threshold = 0.1;
category_types = { '1', '2', '3' };

%%

% try  
  
[curr_ib_detects, curr_could_fix, curr_did_fix] = process_clip_table( ...
    edf_infos, vid_names, vid_p, bbox_p, screen_dims, conf_threshold ...
  , category_types ...
);

% catch err
%   fprintf( '\n\n\n ||| %d failed', si );
%   throw( err );
% end

per_file_outs{si} = table( curr_ib_detects, curr_could_fix, curr_did_fix ...
  , 'va', {'ib_detects', 'could_fix', 'did_fix'} );

end

file_outs = vertcat( per_file_outs{:} );
ib_detects = file_outs.ib_detects;
could_fix = file_outs.could_fix;
did_fix = file_outs.did_fix;

%%

unique_cats = cell( size(ib_detects) );
for i = 1:numel(ib_detects)
  unique_cats{i} = unique( cellfun(@(x) x.category, ib_detects{i}, 'un', 0) );
end

%%  fixation proportions by category (and later split by clip parameters)

fix_props = zeros( 1, size(did_fix, 2) );
for i = 1:size(did_fix, 2)
  fix_props(i) = sum( did_fix(:, i) ) / sum( could_fix(:, i) );
%   fix_props(i) = sum( did_fix(:, i) );
end

subplot( 1, 2, 1 );
bar( fix_props );
set( gca, 'xticklabels', category_names );
ylabel( 'Proportion fixations' );

subplot( 1, 2, 2 );
bar( sum(could_fix, 1) );
set( gca, 'xticklabels', category_names );

%%



%%  number distinct categories per fixation

num_uniques = cellfun( @numel, unique_cats );
hist( num_uniques );

%%

%%

function [all_ib_detects, all_could_fix, all_did_fix] = process_clip_table(...
  edf_infos, vid_names, vid_p, bbox_p, screen_dims, conf_threshold, category_types)

all_ib_detects = [];
all_could_fix = [];
all_did_fix = [];

for i = 1:numel(edf_infos)
  %%
  fprintf( '\n\t %d of %d', i, numel(edf_infos) );
  
  bboxes = load( fullfile(bbox_p, sprintf('%s-bbox', vid_names(i)), 'all_bboxes.mat') );
  vid_reader = VideoReader( fullfile(vid_p, vid_names(i)) );
  im_dims = [ vid_reader.Width, vid_reader.Height ];
  
  %%
  
  fixs = edf_infos{i}.fixations;
  ib_detects = cell( numel(fixs), 1 );
  
  could_fix = zeros( numel(fixs), numel(category_types) );
  did_fix = zeros( size(could_fix) );
  
  for j = 1:numel(fixs)
    fix = fixs(j);
    edf_px = fix.position(1);
    edf_py = fix.position(2);
    
    un_frames = unique( fix.video_frame(~isnan(fix.video_frame)) );
    
    detects = bboxes.detections(un_frames);
    % remove empties
    detects(cellfun(@(x) isequal(x, struct), detects)) = [];  
    if ( ~isempty(detects) )
      detects = horzcat( detects{:} );
    end
    % remove low confidence
    confs = cellfun( @(x) x.conf, detects );
    detects(confs < conf_threshold) = [];
    
    rects = cellfun( ...
        @(x) bbox_to_pixel_rect(x.bbox, im_dims, screen_dims) ...
      , detects, 'un', 0 );
    is_ib = @(r) edf_px >= r(1) & edf_px < r(3) & edf_py >= r(2) & edf_py < r(4);
    ib_rect = cellfun( is_ib, rects );
    
    for k = 1:numel(detects)
      detects{k}.bbox = rects{k};
    end
    
    ib_detects(j) = { detects(ib_rect) };
    
    for k = 1:numel(category_types)
      matches_category = cellfun( @(x) x.category == category_types{k}, detects );
      could_fix(j, k) = any( matches_category );
      did_fix(j, k) = any( ib_rect(matches_category) );
    end
  end
  
  all_ib_detects = [ all_ib_detects; ib_detects ];
  all_could_fix = [ all_could_fix; could_fix ];
  all_did_fix = [ all_did_fix; did_fix ];
  
end

end