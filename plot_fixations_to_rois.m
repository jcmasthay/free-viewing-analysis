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

file_outs = vertcat( per_file_outs{:} );
% ib_detects = file_outs.ib_detects;
% could_fix = file_outs.could_fix;
% did_fix = file_outs.did_fix;

%%

vid_infos = readtable( '~/Downloads/all_data.xlsx' );
vid_infos.code = string( deblank(vid_infos.code) );
vid_infos.prolific_pid = string( deblank(vid_infos.prolific_pid) );

%%

[I, ratings_tbl] = findeach( vid_infos, {'code'} );
ratings_tbl.rating = cellfun( @(x) nanmean(vid_infos.affil_aggr_slider_value(x)), I );

%%  duration of each category

is_prop = true;

if ( is_prop )
  dur = file_outs.dur_did_fix ./ file_outs.dur_could_fix;
else
  dur = file_outs.dur_did_fix;
end

plt_vec = dur;
mask = file_outs.block_type == "A" | file_outs.block_type == "C";
mask(:) = true;

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

%%  stats for duration

[I, C] = findeach( file_outs, {'timestamp', 'identifier', 'block_type', 'category'} ...
  , find(mask) ...
);
means = cellfun( @(x) nanmean(plt_vec(x)), I );

f = fcat.from( C(:, {'block_type', 'category'}) );
anova_outs = dsp3.anovan( means, f, {}, {'block_type', 'category'} );

%%  duration of each category vs. affil-aggr rating

is_prop = true;

if ( is_prop )
  dur = file_outs.dur_did_fix ./ file_outs.dur_could_fix;
else
  dur = file_outs.dur_did_fix;
end

plt_vec = dur;
mask = file_outs.block_type == "A" | file_outs.block_type == "C";

[I, mu_tbl] = findeach( file_outs, {'identifier', 'block_type', 'category'}, mask );
mu_tbl.summary = cellfun( @(x) nanmean(plt_vec(x)), I );
[~, loc] = ismember( mu_tbl.identifier, ratings_tbl.code );
mu_tbl.ratings(:) = nan;
mu_tbl.ratings(loc ~= 0) = ratings_tbl.rating(loc(loc ~= 0));

mask = ~isnan( mu_tbl.ratings ) & ~isnan( mu_tbl.summary );
[I, id, C] = rowsets( 2, mu_tbl, {'block_type'}, 'category' ...
  , 'to_string', 1, 'mask', mask );
[PI, PL] = plots.nest2( id, I, C );
figure(1); clf;
axs = plots.panels( numel(PI) );
for i = 1:numel(axs)
  axes(axs(i));
  [g, v] = ungroupi( PI{i} );
  gscatter( mu_tbl.ratings(v), mu_tbl.summary(v), mu_tbl.category(v) );
  colors = hsv( numel(PI{i}) );
  xlim( axs(i), [-1, 1] );
  for j = 1:numel(PI{i})
    sx = mu_tbl.ratings(PI{i}{j});
    sy = mu_tbl.summary(PI{i}{j});
    ps = polyfit( sx, sy, 1 );
    xs = get( axs(i), 'xlim' );
    hold( axs(i), 'on' );
    h = line( axs(i), xs, polyval(ps, xs) );
    set( h, 'color', colors(j, :) );
    [r, p] = corr( sx, sy );
    txt = sprintf( 'R = %0.3f, P = %0.3f', r, p );
    h = text( axs(i), 0.6, 1 - 0.1 * j, txt );
    set( h, 'color', colors(j, :) );
  end
  title( axs(i), PL{i, 1} );
  xlabel( axs(i), 'affil vs aggressive rating' );
  ylabel( axs(i), 'proportion of time in social roi' );
  shared_utils.plot.match_ylims( axs(i) );
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
is_weighted = false;
only_nz = true;

if ( is_weighted )
  plt_vec = file_outs.weighted_dur_did_fix;
else
  plt_vec = file_outs.dur_did_fix;
end

plt_vec(plt_vec < 0) = nan;

if ( only_nz )
  plt_vec(plt_vec == 0) = nan;
end

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
mask = mask & file_outs.category ~= 'vehicle';
% mask = mask & file_outs.category == 'animal';

plt_type = 'box';
is_violin_or_box = ismember( plt_type, {'box', 'violin'} );
is_box = contains( plt_type, 'box' );

if ( is_violin_or_box )  
  [I, id, C] = rowsets( 3, file_outs ...
    , {} ...
    , {, 'block_type', 'interactive_agency'} ...
    , {'affiliativeness'}, 'mask', mask, 'to_string', true );
  C = strrep( C, '_', ' ');
  
  fi = findeach( id, 1 );
  all_axs = {};
  figs = cell( numel(fi), 1 );
  titles = cell( size(figs) );
  for i = 1:numel(fi)
    figure(i); clf;
    if ( is_box )
      [PI, PL] = plots.nest2( id(fi{i}, 2:end), I(fi{i}), C(fi{i}, 2:end) );
      axs = plots.panels( numel(PI) );
      plots.simple_boxsets( axs, plt_vec, PI, PL );
    else
      axs = plots.violins( plt_vec, I(fi{i}), id(fi{i}, 2:end), C(fi{i}, 2:end) );    
    end
    all_axs{end+1, 1} = axs;
    figs{i} = gcf;
    titles{i} = char( get(get(axs(1), 'title'), 'string'));
  end
  all_axs = vertcat( all_axs{:} );
  if ( is_prop ), ylim( all_axs, [0, 1] ); end
  if ( ~is_prop ), ylim( all_axs, [-100, 1500] ); end
  if ( is_prop ), ylab = 'Proportion of fixation spent in ROI'; end
  if ( ~is_prop ), ylab = 'Fixation duration (ms) spent in ROI'; end
  ylabel( all_axs, ylab );
  
  if ( 1 )
    for i = 1:numel(figs)
      fname = titles{i};
      if ( only_nz ), fname = sprintf( 'only_nonzero_%s', fname ); end
      shared_utils.plot.save_fig( figs{i} ...
        , fullfile(fv_data_directory, 'plots', plt_type, fname), {'png', 'svg', 'epsc'}, true );
    end
  end
  
else
  [I, id, C] = rowsets( 4, file_outs ...
    , {'block_type', 'category'}, {'affiliativeness'}, {'interactive_agency'}, {} ...
    , 'mask', mask, 'to_string', true );
  C = strrep( C, '_', ' ');

  figure(1); clf;
  [axs, hs, xs] = plots.simplest_barsets( plt_vec, I, id, C ...
    , 'summary_func', @nanmean ...
    , 'error_func', @plotlabeled.nansem ...
    , 'add_points', false ...
  );

  set( axs, 'xticklabelrotation', 10 );
  shared_utils.plot.match_ylims( axs );
  if ( is_prop ), ylim( axs, [0, 1] ); end
  
  if ( is_prop ), ylab = 'Proportion of fixation spent in ROI'; end
  if ( ~is_prop ), ylab = 'Fixation duration (ms) spent in ROI'; end
  ylabel( axs, ylab );
  
  fname = char( get(get(axs(1), 'title'), 'string'));
  if ( only_nz ), fname = sprintf( 'only_nonzero_%s', fname ); end
  shared_utils.plot.save_fig( gcf ...
    , fullfile(fv_data_directory, 'plots', 'bar', fname), {'png'}, true );
end

tot_cats = { 'block_type', 'category', 'affiliativeness', 'interactive_agency' };
f = fcat.from( cellstr(file_outs{:, tot_cats}), tot_cats );
anova_mask = intersect( find(mask), findnone(f, 'neutral') );
anova_mask = find( mask );
anova_stats = dsp3.anovan( ...
  plt_vec, f, {'category', 'block_type'}, {'affiliativeness', 'interactive_agency'} ...
  , 'mask', anova_mask ...
);

%%

poss_affil = { 'aggressive', 'affiliative', 'neutral' };
pairs = { [1, 2], [1, 3], [2, 3] };

tot_tbls = table;

for i = 1:numel(pairs)
  a = poss_affil{pairs{i}(1)};
  b = poss_affil{pairs{i}(2)};  
  rs_stats = dsp3.ranksum( ...
    plt_vec, f, {'block_type', 'interactive_agency'}, a, b ...
    , 'mask', anova_mask );
  miss = cellfun( @(x) isnan(x.p), rs_stats.rs_tables );
  keep_tbls = rs_stats.rs_tables(~miss);
  keep_ls = rs_stats.rs_labels(find(~miss));
  keep_tbls = vertcat( keep_tbls{:} );
  setcat( keep_ls, 'affiliativeness', sprintf('%s vs %s', a, b) );
  keep_tbls.affiliativeness = keep_ls(:, 'affiliativeness');
  keep_tbls.interaction = keep_ls(:, 'interactive_agency');
  tot_tbls = [ tot_tbls; keep_tbls ];
end

% rs_stats = dsp3.ranksum( ...
%   plt_vec, f, {'block_type', 'interactive_agency'}, 'aggressive', 'neutral' ...
%   , 'mask', anova_mask )

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