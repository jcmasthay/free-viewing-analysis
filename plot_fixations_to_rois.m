%%  load in data

root_p = fv_data_directory;
fix_info_files = shared_utils.io.findmat( fullfile(root_p, 'fix_infos'));
ct_p = fullfile( root_p, 'edf_samples' );
ct_files = shared_utils.io.findmat( ct_p );
ct_names = shared_utils.io.filenames( ct_files, true );

keep_files = ismember( shared_utils.io.filenames(fix_info_files, true), ct_names );
fix_info_files = fix_info_files(keep_files);

category_names = [ "animal", "human", "vehicle" ];
per_file_outs = cell( numel(fix_info_files), 1 );

for si = 1:numel(fix_info_files)
  fprintf( '\n %d of %d', si, numel(fix_info_files) );
  
  file_outs = shared_utils.io.fload( fix_info_files{si} );
  if ( isempty(file_outs) ), continue; end
  
  fname = shared_utils.io.filenames( fix_info_files{si}, true );
  sesh_dt = datetime( strrep(fname(1:20), '_', ':') );
  sesh_sesh = datestr( sesh_dt, 'mmddyyyy' );
  
  clip_table = shared_utils.io.fload( fullfile(ct_p, fname) );
  clip_table = convert_char_vars_to_string( clip_table );
  file_outs = [ file_outs, clip_table(file_outs.src_index, :) ];
  file_outs.file_name = repmat( string(fname), rows(file_outs), 1 );
  file_outs.session = repmat( string(sesh_sesh), rows(file_outs), 1 );
  
  ind_sets = repmat( 1:rows(file_outs), 1, numel(category_names) );
  ind_ind = repelem( 1:numel(category_names), rows(file_outs) );

  per_cat_outs = file_outs(ind_sets, :);
  per_cat_outs.category = columnize( category_names(ind_ind) );
  
  flat_vars = { ...
    'could_fix', 'did_fix', 'dur_could_fix', 'dur_did_fix', 'weighted_dur_did_fix', 'area_weight' };
  
  for j = 1:numel(flat_vars)
    per_cat_outs.(flat_vars{j}) = columnize( file_outs{:, flat_vars(j)} );
  end
  
  
%   per_cat_outs.dur_could_fix = file_outs.dur_could_fix(:);
%   per_cat_outs.dur_did_fix = file_outs.dur_did_fix(:);
  
  per_file_outs{si} = per_cat_outs;
end

%%

vid_infos = readtable( '~/Downloads/all_data.xlsx' );
vid_infos.code = string( deblank(vid_infos.code) );
vid_infos.prolific_pid = string( deblank(vid_infos.prolific_pid) );

%%

file_outs = vertcat( per_file_outs{:} );
% ib_detects = file_outs.ib_detects;
% could_fix = file_outs.could_fix;
% did_fix = file_outs.did_fix;

%%  duration of each category

is_prop = true;

if ( is_prop )
  dur = file_outs.dur_did_fix ./ file_outs.dur_could_fix;
else
  dur = file_outs.dur_did_fix;
end

plt_vec = dur;
mask = file_outs.block_type == "A" | file_outs.block_type == "C";

[I, id, C] = rowsets( 4, file_outs ...
  , {'block_type'}, {}, {'category'}, {} ...
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

if ( is_prop )
  ylabel( axs(1), 'delta (proportion of fixation)' );
else
  ylabel( axs(1), 'delta (fixation duration)' );
end


%%  duration of each category, affil vs aggressive

is_prop = false;

if ( is_prop )
  dur = file_outs.dur_did_fix ./ file_outs.dur_could_fix;
else
  dur = file_outs.dur_did_fix;
end

sub_each = {'block_type', 'interactive_agency', 'category', 'session'};
[I, C] = findeach( file_outs, sub_each );
C.affiliativeness = repmat( "affil - aggressive", size(I) );

affil_ind = find( file_outs.affiliativeness == 'affiliative' );
aggr_ind = find( file_outs.affiliativeness == 'aggressive' );

mu_diff = nan( size(I) );

for i = 1:numel(I)
  affil = intersect( affil_ind, I{i} );
  aggr = intersect( aggr_ind, I{i} );
  mu_diff(i) = nanmean( dur(affil) ) - nanmean( dur(aggr) );
end

plt_vec = mu_diff;
mask = C.block_type == 'A';

[I, id, C] = rowsets( 4, C ...
  , {'block_type', 'affiliativeness'}, {'interactive_agency'}, {'category'}, {} ...
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

if ( is_prop )
  ylabel( axs(1), 'delta (proportion of fixation)' );
else
  ylabel( axs(1), 'delta (fixation duration)' );
end

%%  duration of each category per fixation

is_prop = false;
is_weighted = true;

if ( is_weighted )
  plt_vec = file_outs.weighted_dur_did_fix;
else
  plt_vec = file_outs.dur_did_fix;
end

plt_vec(plt_vec < 0) = nan;

if ( is_prop )
  if ( is_weighted )
    plt_vec = plt_vec ./ file_outs.area_weight;
  else
    plt_vec = plt_vec ./ file_outs.dur_could_fix;
  end
  plt_vec(~isfinite(plt_vec)) = 0;
end

mask = rowmask( plt_vec );
mask = file_outs.block_type == 'A';

[I, id, C] = rowsets( 4, file_outs ...
  , {'block_type', 'affiliativeness'}, {'interactive_agency'}, {'category'}, {} ...
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
if ( is_prop ), ylim( axs, [0, 1] ); end

%%  clip level proportions of vehicles / humans / animal

%%  within clip, durations of each category

%%  latency to switch to new category
%
%   onset / offset of each category
%     latency to enter category
%
%%  latency to look at *each* roi category, given new onset of *some* roi
%
%
%   

%%

prop_fixated = file_outs.all_dur_did_fix ./ file_outs.all_dur_could_fix;
prop_fixated(isnan(prop_fixated)) = 0;

bar( mean(prop_fixated) );
set( gca, 'xticklabels', category_names );

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
ylim( [0, 1] );

subplot( 1, 2, 2 );
bar( sum(could_fix, 1) );
set( gca, 'xticklabels', category_names );
title( gca, 'Number of fixations that had a category that could be fixated on' );

%%  number distinct categories per fixation

num_uniques = cellfun( @numel, unique_cats );
hist( num_uniques );