load( '~/Downloads/data.mat' );

%%

load( '/Volumes/external3/data/changlab/jamie/free-viewing/data/all_traces.mat' );

%%

samps_p = fullfile( fv_data_directory, 'edf_samples' );
bbox_p = fullfile( fv_data_directory, 'detections' );
samps_files = shared_utils.io.findmat( samps_p );

%%

%%  load survey data

vid_infos = readtable( '~/Downloads/all_data.xlsx' );
vid_infos.code = string( deblank(vid_infos.code) );
vid_infos.prolific_pid = string( deblank(vid_infos.prolific_pid) );
% vid_infos.timestamp = datetime( vid_infos.timestamp );

%%  load video info

video_info = shared_utils.io.fload( ...
  '/Volumes/external3/data/changlab/jamie/free-viewing/videos/vid_info.mat' );

%%

lum_table = shared_utils.io.fload( '/Volumes/external3/data/changlab/jamie/free-viewing/videos/lum.mat' );

%%  load all files

samps = cell( numel(samps_files), 1 );
for i = 1:numel(samps_files)
  fprintf( '\n %d of %d', i, numel(samps_files) );
  samps{i} = shared_utils.io.fload( samps_files{i} );
end
samps = cate1( samps );

%%  find shot times from survey data

samps.search_start_stop = get_merged_shot_starts_stops( samps );
shots = unique( samps(:, {'identifier', 'start', 'stop', }) );

%%  compute traces (e.g. pupil size) over course of movie watching, along with survey data

[I, C] = findeach( samps, {'timestamp', 'Code', 'block_type'} );

of_interest = intersect(unique(C.Code), unique(vid_infos.code));
of_interest = intersect(of_interest, lum_table.code);
% of_interest = {'Odd_Bird_Seduction_Techniques'};

sel = find( ismember( C.Code, of_interest ) );

keep_vars = { 'identifier', 'block_type', 'Code', 'timestamp' };

all_traces = table();
for i = 1:numel(sel)
  fprintf( '\n %d of %d', i, numel(sel) );
  traces = compute_traces( samps(I{sel(i)}, :), vid_infos, keep_vars );
  rest_vars = setdiff( traces.Properties.VariableNames, keep_vars );
  rest_traces = cellfun( @(x) {traces.(x)'}, rest_vars );
  kept_vars = unique( traces(:, keep_vars) );
  assert( rows(kept_vars) == 1 );
  for j = 1:numel(rest_vars)
    kept_vars.(rest_vars{j}) = rest_traces(j);
  end
  all_traces = [ all_traces; kept_vars ];
end

all_traces = sort_traces(all_traces);

%%  time in vs. out of social roi

%{

0. show time spent in social ROI vs. nonsocial ROI
1. compute continuous traces of in vs. out of social ROI
2. select shots with a ~similar baseline probability of being in a social
  ROI (determined by area of rois on screen)
3. 

%}

for i = 1:size(all_traces, 1)
  
    
end

%%  decompose clip-level traces into shot-level traces

%{

1. align pupil traces to each shot
2. filter for shots at least N seconds long
3. compute slope or change in pupil, say, over last second of shot
4a. correlate pupil changes with affil/aggr ratings for the shots
4b. correlate pupil changes with changes in affil/aggr ratings for the shots

%}

shot_level_traces = table();
for i = 1:size(mean_traces, 1)
  fprintf( '\n %d of %d', i, size(mean_traces, 1) );
  
  pup_trace = mean_traces.mean_trace_minus_lum{i};
  rating_trace = mean_traces.rating{i};  
  
  time = mean_traces.time{i};
  code = mean_traces.code(i);
  match_shot = find( ismember(shots.identifier, code) );
  
  get_t = @(t) min( abs(time - t) );
  
  for j = 1:numel(match_shot)
    % for each shot ...
    t0 = shots.start(match_shot(j));
    t1 = shots.stop(match_shot(j));
    
    [~, ind0] = get_t( t0 );
    [~, ind1] = get_t( t1 );
    
    if ( j > 1 )
      [~, prev_i0] = get_t( shots.start(match_shot(j-1)) );
      [~, prev_i1] = get_t( shots.stop(match_shot(j-1)) );
      prev_shot_rating = median( rating_trace(prev_i0:prev_i1) );
    else
      prev_shot_rating = nan;
    end
    
    subset_trace = pup_trace(ind0:ind1);
    shot_rating = median( rating_trace(ind0:ind1) );
    diff_shot_rating = shot_rating - prev_shot_rating;
    was_sign_change = sign( shot_rating ) ~= sign( prev_shot_rating );
    
    shot_level_trace = [table( ...
      {subset_trace}, shot_rating, diff_shot_rating, was_sign_change ...
      , 'va', {'pupil_trace', 'shot_rating', 'change_in_shot_rating', 'was_sign_change'}) ...
      , mean_traces(i, :) ...
      , shots(match_shot(j), {'start', 'stop'}) ];
    shot_level_traces = [ shot_level_traces; shot_level_trace ];
  end
end

%%

vb = 'change_in_shot_rating';
% vb = 'shot_rating';
is_abs_rating = false;
req_sign_change = true;
only_aggr_rating = false;
only_affil_rating = true;

% shot_durs = shot_level_traces.stop - shot_level_traces.start;
num_samps = cellfun( @numel, shot_level_traces.pupil_trace );
ok_shots = num_samps >= 3000;
% ok_shots = ok_shots & num_samps <= 5000;
ok_shots = find( ok_shots );

% compute slope or change in pupil, say, over last second of shot
off0 = 1500;
off1 = off0 + 1500;

shot_subsets = cellfun( @(x) x(off0:off1) ...
  , shot_level_traces.pupil_trace(ok_shots), 'un', 0 );

was_sign_change = shot_level_traces.was_sign_change(ok_shots);

shot_pupil_gradients = cellfun( @(x) mean(gradient(x)), shot_subsets );
% shot_pupil_gradients = cellfun( @(x) x(end) - x(1), shot_subsets );
% shot_ratings = shot_level_traces.shot_rating(ok_shots);
shot_ratings = shot_level_traces.(vb)(ok_shots);

% remove pupil gradients > 5 sds of mean (+/-)
[mu, sigma] = deal( mean(shot_pupil_gradients), std(shot_pupil_gradients) );
lb_thresh = mu - sigma * 5;
ub_thresh = mu + sigma * 5;

surviving_shots = true( numel(shot_pupil_gradients), 1 );

surviving_shots = surviving_shots & strcmp( shot_level_traces.block(ok_shots), 'A' );

surviving_shots = surviving_shots & ...
  shot_pupil_gradients < ub_thresh & ...
  shot_pupil_gradients > lb_thresh;

if ( only_aggr_rating )
  surviving_shots = surviving_shots & shot_ratings < 0;
end
if ( only_affil_rating )
  surviving_shots = surviving_shots & shot_ratings > 0;
end

if ( req_sign_change )
  surviving_shots = surviving_shots & was_sign_change;
end

surviving_shots = surviving_shots & ~isnan( shot_ratings );

[shot_ratings, shot_pupil_gradients] = deal( ...
  shot_ratings(surviving_shots), shot_pupil_gradients(surviving_shots) );

if ( is_abs_rating )
  shot_ratings = abs( shot_ratings );
end

if ( 1 )
  [I, C] = findeach( ...
    shot_level_traces(ok_shots(surviving_shots), :), {} );
  C = {C};
else
  [I, C] = findeach( ...
    shot_level_traces(ok_shots(surviving_shots), :), 'code' );
end

figure(1); clf;

axs = plots.panels( numel(I) );

for i = 1:numel(I)
  
axes( axs(i) );
ind = I{i};

% plot( shot_pupil_gradients{2} );
scatter( shot_ratings(ind), shot_pupil_gradients(ind) );

if ( i == 1 )
  xlabel( strrep(vb, '_', ' ') );
  ylabel( sprintf(...
    'Mean pupil gradient from [%d, %d]ms relative to shot', off0, off1) );
end

t_str = strjoin( C{i, :} );
[r, p] = corr( shot_ratings(ind), shot_pupil_gradients(ind) );
title( compose("%s \n | R = %0.3f, P = %0.3f", t_str, r, p) );

ps = polyfit( shot_ratings(ind), shot_pupil_gradients(ind), 1 );
lims = get( gca, 'xlim' );
ys = polyval( ps, lims );
line( lims, ys );
% xlim( [-1, 1] );
ylim( [-2e-3, 2e-3] );

end

%%
%   resample traces to common time axis
resamp_I = findeach( all_traces, 'identifier' );
for i = 1:numel(resamp_I)
    fprintf( '\n %d of %d', i, numel(resamp_I) );
    ri = resamp_I{i};

    max_t = min( arrayfun(@(x) max(all_traces.video_t{x}), ri) );
    min_t = max( arrayfun(@(x) min(all_traces.video_t{x}), ri) );
    round_to = 1e3;
    time_axis = min_t:1/round_to:(floor(max_t * round_to) / round_to);

    aligned_traces = nan( numel(ri), numel(time_axis) );
    aligned_lums = nan( numel(ri), numel(time_axis) );
    aligned_ratings = nan( numel(ri), numel(time_axis) );
    aligned_pos = nan( numel(ri), numel(time_axis), 2 );

    clip_name = all_traces.Code{ri(1)};
    if contains(clip_name, 'Mcq')
        if contains(clip_name, 'Huddle') | ...
            contains(clip_name, 'Swimmin')
            fps = 30;
        elseif contains(clip_name, 'Feed') | ...
                contains(clip_name, 'Pool')
            fps = 25;
        else
            fps = 24;
        end
    else
        fps = 30;
    end

    for j = 1:numel(ri)
        fprintf( '\n\t %d of %d', j, numel(ri) );
        video_t = all_traces.video_t{ri(j)};
        pupil_size = all_traces.cleaned_pup{ri(j)};
        % pupil_size = all_traces.smoothed_pup{ri(j)};
        ratings = all_traces.affil_aggr_ratings{ri(j)};
        pos = all_traces.position{ri(j)};

        lum = lum_table.luminance{string(lum_table.code) == all_traces.Code(ri(j))};
        lum_t = (0:numel(lum)-1)/fps;

        rounded_video_t = floor( video_t * round_to ) / round_to;

        for k = 1:numel(rounded_video_t)
            [~, ind] = min( abs(rounded_video_t(k) - time_axis) );
            aligned_traces(j, ind) = pupil_size(k); 
            aligned_ratings(j, ind) = ratings(k);
            aligned_pos(j, ind, :) = pos(:, k);
        end
        
        if j > 1
            aligned_lums(j, :) = aligned_lums(1, :);
        else
            for k = 1:numel(lum_t)
                [~, ind] = min( abs(lum_t(k) - time_axis) );
                aligned_lums(j, ind) = lum(k);
            end
        end

        aligned_traces(j, :) = fillmissing( aligned_traces(j, :), 'linear' );
        aligned_lums(j, :) = fillmissing( aligned_lums(j, :), 'linear' );
        aligned_ratings(j, :) = fillmissing( aligned_ratings(j, :), 'linear' );
        for k = 1:size(aligned_pos, 3)
          aligned_pos(j, :, k) = fillmissing( aligned_pos(j, :, k), 'linear' );
        end
    end

    smoothed_lums = smoothdata( aligned_lums, 2, 'gaussian', 5e3 );
    smoothed_traces = smoothdata( aligned_traces, 2, 'gaussian', 5e3 );
    
    for p = 1:numel( ri )
        all_traces.aligned_lum( ri(p) ) = {smoothed_lums( p, : )};
        all_traces.aligned_trace( ri(p) ) = {smoothed_traces( p, : )};
        all_traces.aligned_rating( ri(p) ) = {aligned_ratings( p, : )};
        all_traces.aligned_time( ri(p) ) = {time_axis};
        all_traces.aligned_position( ri(p) ) = {squeeze(aligned_pos(p, :, :))};
    end
end

%%

all_traces.within_screen_proportion = nan( size(all_traces, 1), 1 );

for i = 1:size(all_traces, 1)

position_trace0 = all_traces.aligned_position{i};
  
[ibx, iby] = deal(...
    position_trace0(:, 1) >= 0 & position_trace0(:, 1) <= 1600 ...
  , position_trace0(:, 2) >= 0 & position_trace0(:, 2) <= 900 ...
);

ib = ibx & iby;

all_traces.within_screen_proportion(i) = pnz( ib );

end

m = ~strcmp( all_traces.block_type, 'C' );
bins = linspace( 0, 1, 20 );
figure(1); clf; hist(all_traces.within_screen_proportion(m),bins);

% median( all_traces.within_screen_proportion(m) )
% pnz( all_traces.within_screen_proportion(m) > 0.9 )

%%  remove first N ms of position data of each shot by filling them with NaNs

all_traces.aligned_truncated_position = all_traces.aligned_position;

remove_first_ms = 500;

for i = 1:size(all_traces, 1)

fprintf( '\n %d of %d', i, size(all_traces, 1) );

[boundaries, shot_indices] = find_shot_boundaries( ...
  shots, all_traces.aligned_time{i}, all_traces.Code{i} );

for j = 1:size(boundaries, 1)
  [i0, i1] = deal( boundaries(j, 1), boundaries(j, 2) );
  i1 = min( i0 + remove_first_ms, i1 );
  all_traces.aligned_truncated_position{i}(i0:i1, :) = nan;
end

end

%%



%%

% clip_id = "LongTailedMacaque_Human_FeedToolUsePlay";
% clip_id = "Orangutan_Human_TrainforWild";
clip_id = "Marmoset_Harvest_Groom";

[I, C] = findeach( corr_tbls, {'identifier', 'block_type'} ...
  , strcmp(corr_tbls.block_type, 'A') )
C.r = cellfun( @(x) mean(corr_tbls.r(x)), I );
sortrows( C, 'r' )

twin = [15, 30];

search_for = find( ...
  all_traces.identifier == clip_id & ...
  strcmp(all_traces.block_type, 'A') );

% search_for = [ search_for_a(1), search_for_b(1) ];
% search_for = search_for_b;

figure(1); clf;
% axs = plots.panels( numel(search_for) );

axs = plots.panels( 2 );

for i = 1:numel(search_for)
  
ind_x = 1;
ind_y = 2;

time_axis = all_traces.aligned_time{search_for(i)};
pos_trace = all_traces.aligned_position{search_for(i)};

win = 250;
time_axis = time_axis(1:win:end);
pos_trace = pos_trace(1:win:end, :);

plot( axs(ind_x), time_axis, pos_trace(:, 1), 'DisplayName', sprintf('viewing x %d', i) );
plot( axs(ind_y), time_axis, pos_trace(:, 2), 'DisplayName', sprintf('viewing y %d', i) );

hold( axs(ind_x), 'on' );
hold( axs(ind_y), 'on' );
legend;

title( axs(ind_x), plots.strip_underscore(compose("%s x", clip_id)) );
title( axs(ind_y), plots.strip_underscore(compose("%s y", clip_id)) );

end

shared_utils.plot.match_ylims( axs );
xlim( axs, twin );

%%  compute correlations between scan paths on multiple viewings

mask = all_traces.within_screen_proportion >= 0.85;
% mask(:) = true;

downsample_step = 50;

% pos_var = 'aligned_position';
pos_var = 'aligned_truncated_position';

[I, C] = findeach( all_traces, {'identifier', 'block_type'}, mask );

corr_tbls = table();
for i = 1:numel(I)
  fprintf( '\n %d of %d', i, numel(I) );
  
  inds = I{i};
  pairs = bfw.pair_combination_indices( numel(inds) );
  
  for j = 1:size(pairs, 1)
    [i0, i1] = deal( pairs(j, 1), pairs(j, 2) );
    ind0 = inds(i0);
    ind1 = inds(i1);
    
    position_trace0 = all_traces.(pos_var){ind0};
    position_trace1 = all_traces.(pos_var){ind1};
    
    if ( 1 )
      position_trace0 = position_trace0(1:downsample_step:end, :);
      position_trace1 = position_trace1(1:downsample_step:end, :);
    end
    
    [ibx0, iby0] = deal(...
        position_trace0(:, 1) >= 0 & position_trace0(:, 1) <= 1600 ...
      , position_trace0(:, 2) >= 0 & position_trace0(:, 2) <= 900 ...
    );
  
    [ibx1, iby1] = deal(...
        position_trace1(:, 1) >= 0 & position_trace1(:, 1) <= 1600 ...
      , position_trace1(:, 2) >= 0 & position_trace1(:, 2) <= 900 ...
    );
  
    keep_pos = ibx0 & iby0 & ibx1 & iby1;
    
    [trace0_x, trace0_y] = deal(...
        zscore(position_trace0(keep_pos, 1)) ...
      , zscore(position_trace0(keep_pos, 2)) ...
    );
  
    position_trace0 = [ trace0_x; trace0_y ];
  
    [trace1_x, trace1_y] = deal(...
        zscore(position_trace1(keep_pos, 1)) ...
      , zscore(position_trace1(keep_pos, 2)) ...
    );
  
    position_trace1 = [ trace1_x; trace1_y ];
    
    time_trace0 = all_traces.aligned_time{ind0}(1:downsample_step:end);
    time_trace1 = all_traces.aligned_time{ind1}(1:downsample_step:end);
    
    time_trace0 = time_trace0(keep_pos);
    time_trace1 = time_trace1(keep_pos);
    
    [r, p] = corr( position_trace0, position_trace1, 'rows', 'complete' );
    corr_tbl = [ C(i, :), table(r, p, i0, i1, 'VariableNames', {'r', 'p', 'i0', 'i1'}) ];
    corr_tbls = [ corr_tbls; corr_tbl ];
    
    if ( 0 )
      figure(1); clf;
      ax = subplot( 2, 1, 1 );
      hold( ax, 'on' );
      plot( time_trace0, position_trace0(1:numel(time_trace0/2)) );
      plot( time_trace1, position_trace1(1:numel(time_trace1/2)) );
      ax = subplot( 2, 1, 2 );
      hold( ax, 'on' );
      plot( time_trace0, position_trace0(numel(time_trace0/2)+1:end) );
      plot( time_trace1, position_trace1(numel(time_trace1/2)+1:end) );

      title( ax, plots.strip_underscore(sprintf('%s_%s', C.identifier{i}, C.block_type{i})) );
      save_p = fullfile( fv_data_directory, 'plots/traces/scanpath_pairs' );
      shared_utils.io.require_dir( save_p );
      saveas( gcf, fullfile(save_p ...
        , sprintf('%s_%s_%d_%d.fig', C.identifier{i}, C.block_type{i}, i0, i1)) );
    end
  end
end

%%

[I, C] = findeach( corr_tbls, 'identifier' );

C.ps = nan( numel(I), 1 );
for i = 1:numel(I)
  ind_a = intersect( I{i}, find(strcmp(corr_tbls.block_type, 'A')) );
  ind_b = intersect( I{i}, find(strcmp(corr_tbls.block_type, 'B')) );
  
  if ( isempty(ind_a) || isempty(ind_b) )
    continue
  end
  
  C.ps(i) = ranksum( corr_tbls.r(ind_a), corr_tbls.r(ind_b) );
end

C(C.ps < 0.05, :)

%%  compare scan path correlation coefficients between block types

[I, C] = findeach( corr_tbls, {'identifier', 'i0', 'i1'} );
a_ind = find( strcmp(corr_tbls.block_type, 'A') );
b_ind = find( strcmp(corr_tbls.block_type, 'B') );
c_ind = find( strcmp(corr_tbls.block_type, 'C') );

coeff_pairs = table();
for i = 1:numel(I)
  ind = I{i};
  curr_a_ind = intersect( a_ind, ind );
  curr_b_ind = intersect( b_ind, ind );
  curr_c_ind = intersect( c_ind, ind );
  
  if ( isempty(curr_a_ind) ), ra = nan; else; ra = corr_tbls.r(curr_a_ind); end
  if ( isempty(curr_b_ind) ), rb = nan; else; rb = corr_tbls.r(curr_b_ind); end
  if ( isempty(curr_c_ind) ), rc = nan; else; rc = corr_tbls.r(curr_c_ind); end
  
  pair_a_b = table( ra, rb, {'A'}, {'B'} ...
    , 'va', {'set1', 'set2', 'block_type1', 'block_type2'} );
  pair_a_c = table( ra, rc, {'A'}, {'C'} ...
    , 'va', {'set1', 'set2', 'block_type1', 'block_type2'} );
  pair_b_c = table( rb, rc, {'B'}, {'C'} ...
    , 'va', {'set1', 'set2', 'block_type1', 'block_type2'} );
  
  pair_a_b = [ pair_a_b, C(i, :) ];
  pair_a_c = [ pair_a_c, C(i, :) ];
  pair_b_c = [ pair_b_c, C(i, :) ];
  
  coeff_pairs = [ coeff_pairs; pair_a_b; pair_a_c; pair_b_c ];
end

[I, coeff_pair_means] = findeach( ...
  coeff_pairs, {'block_type1', 'block_type2', 'identifier'} );
[coeff_pair_means.set1, coeff_pair_means.set2] = cellfun(...
  @(x) deal(nanmean(coeff_pairs.set1(x)), nanmean(coeff_pairs.set2(x))), I);

[I, coeff_pair_stats] = findeach( coeff_pair_means, {'block_type1', 'block_type2'} );

coeff_pair_stats.ps = nan( numel(I), 1 );
coeff_pair_stats.stats = cell( numel(I), 1 );
for i = 1:numel(I)
  try
    ind = I{i};
%     [coeff_pair_stats.ps(i), ~, coeff_pair_stats.stats{i}]
    [~, coeff_pair_stats.ps(i)] = ttest( ...
      coeff_pair_means.set1(ind), coeff_pair_means.set2(ind) );
  end
end

[I, clip_idents] = findeach( coeff_pairs, {'identifier'} ...
  , strcmp(coeff_pairs.block_type1, 'A') & ...
    strcmp(coeff_pairs.block_type2, 'B') );
clip_idents.r1 = cellfun( @(x) mean(coeff_pairs.set1(x)), I );
clip_idents.r2 = cellfun( @(x) mean(coeff_pairs.set2(x)), I );
sortrows( clip_idents, 'r1', 'descend' )

%%  compute all pair-wise correlations between scan paths on different viewings

% pos_var = 'aligned_position';
pos_var = 'aligned_truncated_position';

[I, C] = findeach( all_traces, {'block_type'} );
clip_lens = cellfun( @numel, all_traces.aligned_time );
min_len = min( clip_lens );

trace_corrs = table();
for i = 1:numel(I)
  fprintf( '\n %d of %d', i, numel(I) );
  ind = I{i};
  mat_r = nan( numel(ind) );
  [~, ord] = sort( all_traces.identifier(ind) );
  ind = ind(ord);
  for j = 1:numel(ind)
    fprintf( '\n\t %d of %d', j, numel(ind) );
    for k = 1:numel(ind)
      ind0 = ind(j);
      ind1 = ind(k);
      position_trace0 = all_traces.(pos_var){ind0};
      position_trace1 = all_traces.(pos_var){ind1};
      position_trace0 = columnize( position_trace0(1:min_len, :) );
      position_trace1 = columnize( position_trace1(1:min_len, :) );
      
      [r, p] = corr( position_trace0, position_trace1 );
      mat_r(j, k) = r;      
%       trace_corrs(end+1, :) = [ C(i, :), table(r, p, j, k, 'va', {'r', 'p', 'j', 'k'}) ];
    end
  end
  
  idents = all_traces.identifier(ind);
  trace_corrs(end+1, :) = [ C(i, :), table({mat_r}, {idents}, 'va', {'r', 'identifiers'}) ];
end

%%

[I, C] = findeach( trace_corrs, {'block_type'} );
figure(1); clf;
axs = plots.panels( numel(I) );
for i = 1:numel(I)
  idents = trace_corrs.identifiers{I{i}};
  ind = ismember(...
      idents ...
    , correspondence_tbl.bot_clips{...
    strcmp(correspondence_tbl.block_type, C.block_type{i})});
  
  imagesc( axs(i), trace_corrs.r{I{i}}(ind, ind) );  
  [~, ~, ic] = unique( trace_corrs.identifiers{I{i}}(ind) );
  diffs = find( diff(ic) > 0 );
  
  title( axs(i), C.block_type{i} );
  colorbar( axs(i) );
  shared_utils.plot.set_clims( axs(i), [-0.2, 0.3] );
  hold( axs(i), 'on' );
%   for j = 1:numel(diffs)
%     text( axs(i), 0, diffs(j), trace_corrs.identifiers{I{i}}(diffs(j)) );
%   end
end

%%

target_pair = find( ...
  all_traces.identifier == "Mcq_Mongoose_Play" & ...
  strcmp(all_traces.block_type, 'A') );

tp = target_pair(1:2);

figure(1); clf; hold on;
plot( all_traces.aligned_time{tp(1)}, all_traces.aligned_position{tp(1)}(:, 1) );
plot( all_traces.aligned_time{tp(2)}, all_traces.aligned_position{tp(2)}(:, 1) );

%%  average correlations over pairs of viewings of individual clips

[I, mu_corr_tbls] = findeach( corr_tbls, {'identifier', 'block_type'} );
mu_corr_tbls.r = cellfun( @(x) mean(corr_tbls.r(x)), I );

%%  examine clips with highest correspondence (r value)

[I, C] = findeach( corr_tbls, {'block_type', 'identifier'} );
C.r = cellfun( @(x) nanmean(corr_tbls.r(x)), I );

[I, correspondence_tbl] = findeach( C, 'block_type' );
correspondence_tbl.top_clips(:) = { [] };

for i = 1:numel(I)
  ind = I{i};
  rs = C.r(ind);
  quants = prctile( rs, [10, 90, 100] );
  top_ind = ind(rs >= quants(2));
  bot_ind = ind(rs < quants(1));
  top_clips = C.identifier(top_ind);
  bot_clips = C.identifier(bot_ind);
  correspondence_tbl.top_clips{i} = [top_clips, C.r(top_ind) ];
  correspondence_tbl.bot_clips{i} = [bot_clips, C.r(bot_ind) ];
end

%%

[I, C] = findeach( mu_corr_tbls, {'block_type'} );
figure(1); clf;
axs = plots.panels( numel(I) );
for i = 1:numel(I)
  bins = linspace( -1, 1, 50 );
  hist( axs(i), mu_corr_tbls.r(I{i}), bins );
  xlim( axs(i), [-1, 1] );
  title( axs(i), C.block_type{i} );
  med = median( mu_corr_tbls.r(I{i}) );
  hold( axs(i), 'on' );
  shared_utils.plot.add_vertical_lines( axs(i), med );
end

%%

save_path = '~/Downloads/traces_data.mat';
save(save_path, 'all_traces', '-v7.3');

%%
all_traces = sortrows(all_traces, {'Code', 'block_type'});

[groups, codes, blocktypes] = findgroups(all_traces.Code, all_traces.block_type);

avgPupilSize = splitapply(@(x) {nanmean(cell2mat(x), 1)}, all_traces.aligned_trace, groups);
firstLumVals = splitapply(@(x) x(1), all_traces.aligned_lum, groups);
firstRatings = splitapply(@(x) x(1), all_traces.aligned_rating, groups);
firstTime = splitapply(@(x) x(1), all_traces.aligned_time, groups);

corrected_pup = cellfun( @(x, y) zscore(x) + zscore(y), avgPupilSize, firstLumVals, 'un', 0 );

mean_traces = table(codes, blocktypes, avgPupilSize, corrected_pup, ...
    firstLumVals, firstRatings, firstTime, 'VariableNames',... 
    {'code', 'block', 'mean_trace', 'mean_trace_minus_lum', 'lum', 'rating', 'time'});

%%
codes = unique(mean_traces.code);
mean_diffs = nan(length(codes), 1);
mean_ratings = nan(length(codes), 1);
corr_coefs = nan(length(codes), 1);
corr_pvals = nan(length(codes), 1);

for i = 1:numel(codes)
    code = codes{i};

    a_id = find( ...
      strcmp(mean_traces.code, code) & ...
      strcmp(mean_traces.block, 'A'));
    b_id = find( ...
      strcmp(mean_traces.code, code) & ...
      strcmp(mean_traces.block, 'B'));
    c_id = find( ...
      strcmp(mean_traces.code, code) & ...
      strcmp(mean_traces.block, 'C'));
    
    a_trace = mean_traces.mean_trace{ a_id };
    b_trace = mean_traces.mean_trace{ b_id };
    c_trace = mean_traces.mean_trace{ c_id };
    rating = mean_traces.rating{ a_id };
    lum = mean_traces.lum{ a_id };
    time = mean_traces.time{ a_id };
    
    plt_title = strrep(code, '_', '-');

    figure(i); clf;

    subplot(2, 1, 1);
    plot(time, zscore(movmean(a_trace, 100)), 'displayname', 'pupil (block A)'); hold on;
    plot(time, -zscore(movmean(lum, 100)), 'displayname', 'luminance'); hold on;
    legend;
    title(plt_title);
    
    pup_minus_lum = zscore(movmean(a_trace, 100)) + zscore(movmean(lum, 100));

    subplot(2, 1, 2);
    plot(time, pup_minus_lum, 'displayname', 'pup - lum'); hold on;
    plot(time, 2*(rating), 'displayname', 'rating');
    legend;
    title(plt_title);

%     figure_position = [100, 100, 800, 500];
%     set(gcf, 'Position', figure_position);

    save_path = strcat('/Users/efedogruoz/Desktop/pupsize_plots/', code, '.png');
    %saveas(gcf, save_path);
    %close(gcf);

    [r, p] = corr(pup_minus_lum(1:100:end)', rating(1:100:end)');
    corr_coefs(i) = r;
    corr_pvals(i) = p;

    diff = a_trace - c_trace;

    mean_diff = mean(diff, 2);
    mean_rating = mean(rating, 2);

    mean_diffs(i) = mean_diff;
    mean_ratings(i) = mean_rating;
end

clip_level_stats = table(codes, mean_diffs, mean_ratings,...
    corr_coefs, corr_pvals, 'VariableNames',... 
    {'code', 'mean_a_c_diff', 'mean_rating', 'corr_a_lum', 'pval'});
  
%%
% calculate arousal by behavior category

for p = 1:height( mean_traces )
    mean_traces.arousal(p) = {zscore( mean_traces.mean_trace{p} ) + zscore( mean_traces.lum{p} )};
end

path = '~/Downloads/behavior_shots/';
files = dir(fullfile(path, '*.csv'));

arousal_levels = cell(numel(files), 2);

for i = 1:length(files)
    file_name = files(i).name;
    behavior = file_name(1:end-4);
    table = readtable([ path '/' file_name ]);
    table.arousal_change(:) = nan;
    
    arousal_levels( i, 1 ) = { behavior };
    vec = cell( height( table ), 1 );

    figure(i); clf;

    prev_code = 0;
    prev_stop_time = 0;

    for p = 1:length(vec)
        
        code = table.code{p};
        start_time = table.start(p);
        stop_time = table.stop(p);

        if start_time < 11
            continue
        end

        if 1
            if strcmp(code, prev_code) & (start_time == prev_stop_time)
                prev_code = code;
                prev_stop_time = stop_time;
                continue
            end
        end
        
        scene_ids = find(strcmp(mean_traces.code, code));

        if isempty( scene_ids )
            continue
        end

        scene_id = scene_ids(1);
        idx = ( mean_traces.time{scene_id} >= start_time ) & ...
            ( mean_traces.time{scene_id} <= stop_time );
        
        arousal = mean_traces.arousal{scene_id}(idx);
        
        if 1
            if numel( arousal ) < 3000
                continue
            end
        end

        rating = mean(mean_traces.rating{scene_id}(idx));

        plot( arousal ); hold on;
        vec(p) = { arousal };
        table.arousal_change( p ) = 10000 * max(gradient(arousal));
        prev_code = code;
        prev_stop_time = stop_time;
    end
    title( behavior );
%     close(gcf);

    arousal_levels( i, 2 ) = {vec};
    
%     x = input
end

%%


%%



%%

%{

1. median pupil over shot

%}

%%
behaviors = arousal_levels(:, 1);

for p = 1:numel(behaviors)
    behavior = behaviors{p};
    traces = arousal_levels(p, 2);
    traces = traces{1};

    slopes_vec = [];
    changes_vec = [];

    for i = 1:numel(traces)
        trace = traces{i};
        if numel( trace ) == 0
            continue
        end
        
        slopes_vec = [slopes_vec mean(gradient(trace))];
        change = trace(2000) - trace(1000);
        changes_vec = [changes_vec change];
    end
    
    arousal_levels(p, 3) = {median(changes_vec)};
    arousal_levels(p, 4) = {std(changes_vec)};
    arousal_levels(p, 5) = {slopes_vec};
end

% order = [1 2 3 12 7 4 13 6 5 8 9 10];
order = [1 2 3 12 7 4 6 5 8 9 10];
ord_arousal_levels = arousal_levels(order, :);

figure; clf;
bar(cell2mat(ord_arousal_levels(:,3)));
hold on;
errorbar(1:rows(ord_arousal_levels), cell2mat(ord_arousal_levels(:,3)), cell2mat(ord_arousal_levels(:,4)), 'k', 'linestyle', 'none');
xticks(1:numel(ord_arousal_levels(:,1)));
xticklabels(ord_arousal_levels(:,1));

%%

for i = 1:height(unique(mean_traces.code))
    figure(i); clf;
    id = (i - 1) * 3;
    plot(mean_traces.mean_trace{id + 1}, 'displayname', 'A'); hold on;
    plot(mean_traces.mean_trace{id + 2}, 'displayname', 'B');
    plot(mean_traces.mean_trace{id + 3}, 'displayname', 'C');
    legend;
    title(mean_traces.code{id + 1});
    close(gcf);
end


%%

smoothed_traces = smoothdata( aligned_traces, 2, 'gaussian', 5e3 );
mean_smoothed_pup = nanmean( smoothed_traces, 1 );

smoothed_traces_c = smoothdata( c_aligned_traces, 2, 'gaussian', 5e3 );
meaned_smoothed_pup_c = nanmean( smoothed_traces_c, 1 );

lum_trace = aligned_lums(1, :);
lum_trace = smoothdata( lum_trace, 'gaussian', 5e3 );

rating_trace = aligned_ratings(1, :);

figure(2); clf;

ax = subplot( 2, 1, 1 );
hold( ax, 'on' );
plot( time_axis, -zscore(lum_trace), 'displayname', 'z-scored inverse luminance' ); hold on;
plot( time_axis, zscore(mean_smoothed_pup), 'displayname', 'z-scored pup' );
plot( time_axis, zscore(meaned_smoothed_pup_c), 'displayname', 'z-scored pup c' );
legend;

ax = subplot( 2, 1, 2 );
hold( ax, 'on' );
%diff_trace = zscore( mean_smoothed_pup - meaned_smoothed_pup_c );
diff_trace = zscore( mean_smoothed_pup ) + zscore( lum_trace );
plot( time_axis, diff_trace, 'displayname', 'z-scored pup - pup c' );
plot( time_axis, (rating_trace), 'displayname', 'rating' );
legend;

[r, p] = corr( lum_trace(:), mean_smoothed_pup(:), 'rows', 'complete' );
[r2, p2] = corr( diff_trace(:), rating_trace(:), 'rows', 'complete' );

%%
window_size = 1000;

sm_pup = movmean(mean_smoothed_pup, window_size);
sm_lum = movmean(lum_trace, window_size);
sm_time = movmean(time_axis, window_size);
sm_rating = movmean(rating_trace, window_size/100);

figure(1); clf;
sm_diff = zscore(sm_pup) + zscore(sm_lum);
plot( sm_time, -zscore(sm_lum), 'displayname', 'lum'); hold on;
plot( sm_time, zscore(sm_pup), 'displayname', 'pup');
legend;

figure(2); clf;
plot( sm_time, zscore(sm_diff), 'displayname', 'pup - lum'); hold on;
plot( sm_time, zscore(sm_rating), 'displayname', 'rating');
legend;

corr( sm_diff(:), sm_rating(:), 'rows', 'complete')

%%

figure(1); clf;
axs = plots.panels( [3, 1] );
smooth_ax_ind = 1;
rating_ax_ind = smooth_ax_ind + 1;
lum_ax_ind = rating_ax_ind + 1;

mean_only = false;

title_str = plots.cellstr_join( ...
  {unique(all_traces(:, {'Code'}))} );
title_str = strjoin( title_str, ' | ' );
title_str = plots.strip_underscore( title_str );

ylim( axs(rating_ax_ind), [-1, 1] );

mean_smoothed_pup = nanmean( aligned_traces, 1 );
h = plot( axs(smooth_ax_ind), time_axis, mean_smoothed_pup );
set( h, 'linewidth', 2 );
hold( axs(smooth_ax_ind), 'on' );

for i = 1:size(all_traces, 1)

    if ( ~mean_only ) 
        plot( axs(smooth_ax_ind), all_traces.video_t{i}, all_traces.smoothed_pup{i} );
        hold( axs(smooth_ax_ind), 'on' );
    end
  
  plot( axs(rating_ax_ind), all_traces.video_t{i}, all_traces.affil_aggr_ratings{i} );
  hold( axs(rating_ax_ind), 'on' );
  
  if ( i == 0 )
    plot( axs(rating_ax_ind), all_traces.video_t{i}, all_traces.affil_aggr_ratings{i} - all_traces.affil_aggr_ratings_err{i} );
    plot( axs(rating_ax_ind), all_traces.video_t{i}, all_traces.affil_aggr_ratings{i} + all_traces.affil_aggr_ratings_err{i} );
  end
  
  fps = 30;

  lum = lum_table.luminance{string(lum_table.code) == all_traces.Code(i)};
  smooth_lum1 = smoothdata( lum, 'gaussian', floor(fps * 0.5) );
  smooth_lum2 = smoothdata( lum, 'gaussian', floor(fps * 5) );
 
  hold( axs(lum_ax_ind), 'on' );
  plot( axs(lum_ax_ind), (0:numel(lum)-1)/fps, smooth_lum1(:)' );
  plot( axs(lum_ax_ind), (0:numel(lum)-1)/fps, smooth_lum2(:)' );
  
  title( axs(smooth_ax_ind), title_str );
  title( axs(rating_ax_ind), title_str );
  title( axs(lum_ax_ind), title_str );
  
  ylabel( axs(smooth_ax_ind), 'Pupil size' );
  ylabel( axs(rating_ax_ind), 'Affil vs. aggr rating' );
  ylabel( axs(lum_ax_ind), 'Video luminance' );
  
  shared_utils.plot.match_xlims( axs );
end

offset = all_traces.smoothed_pup{2}(1) - all_traces.smoothed_pup{1}(1);
fprintf( '\n offset: %.1f', offset);

n = 1; % which pupil trace to use for stats

ratings = cell2mat(all_traces.affil_aggr_ratings(1));
pup = cell2mat(all_traces.smoothed_pup(n));
scaling_factor = numel(pup)/numel(lum);
% expanded_lum = interp1(1:numel(lum), lum, linspace(1, numel(lum), scaling_factor*numel(lum)), 'linear');
expanded_lum = interp1(1:numel(lum), smooth_lum2, linspace(1, numel(lum), scaling_factor*numel(lum)), 'linear');

min_rows = min([numel(ratings), numel(pup), numel(expanded_lum)]);

X = pup(1:min_rows)';
Y = ratings(1:min_rows)';
Z = expanded_lum(1:min_rows)';

X = X(1:100:end);
Y = Y(1:100:end);
Z = Z(1:100:end);
t = all_traces.video_t{n}(1:100:end);

if (1)

    disp(title_str);
    fprintf( '\n pup x lum: %.2f', corr(X, Z, 'rows', 'complete') );
    fprintf( '\n pup mean x lum: %.2f', corr(pup_mean, Z, 'rows', 'complete') );
    fprintf( '\n pup x ratings: %.2f', corr(X, Y, 'rows', 'complete') );
    fprintf( '\n lum x ratings: %.2f\n', corr(Y, Z, 'rows', 'complete') );
    
    [r, p] = partialcorr(X, Y, Z, 'rows', 'complete');
    
    fprintf( '\n pup x rating, controlling\n for lum: r = %.2f, p = %0.4f', r, p );

end

%%

X = pup(1:min_rows)';
Y = ratings(1:min_rows)';
Z = expanded_lum(1:min_rows)';

X = X(1:100:end);
Y = Y(1:100:end);
Z = Z(1:100:end);
t = all_traces.video_t{n}(1:100:end);
%t = t(:, 1:end-1);

pupil_size = X;
rating = Y;
luminance = Z;

lum_with_intercept = [luminance(:), ones(length(luminance), 1)];

[beta,bint,residuals] = regress(pupil_size, lum_with_intercept);

pup_pred = luminance*beta(1) + ones(length(luminance), 1) * beta(2);
pup_pred_prime = luminance*beta_prime(1) + ones(length(luminance), 1) * beta_prime(2);

figure(2); clf;

plot( t, pup_pred, 'r', 'displayname', 'predicted pupil size' ); hold on;
plot( t, pupil_size, 'b', 'displayname', 'actual pupil size' ); hold on;
plot( t, pup_pred_prime, 'g', 'displayname', 'univ-predicted pupil size' ); hold on;
legend;

remainder = pupil_size - pup_pred;
remainder_prime = pupil_size - pup_pred_prime;

if (1)
    figure(4); clf;
    
    plot( t, zscore(remainder), 'r', 'displayname', 'remainder' ); hold on;
    %plot( remainder_prime, 'b', 'displayname', 'univ-predicted remainder' ); hold on;
    plot( t, zscore(rating), 'g', 'displayname', 'aggr-affil ratings' )
    legend;
    
    fprintf( '\n r = %.4f', corr(pup_pred, pupil_size) );
end
%%

% Assuming X is the variable you want to regress out, and Y is your dependent variable
% X and Y should be vectors of the same length

pupil_size = zscore( X(1:end-10) );
if ( 0 )
    rating = zscore( abs(Y) );
else
    rating = zscore( Y(1:end-10) );
end
luminance = zscore( Z(1:end-10) );

% Add a column of ones to X to include an intercept in the model
% X_with_intercept = [ones(length(X), 1) x(:)];
X_with_intercept = [pupil_size(:)];

% Regress Y on X to get the residuals
[beta,~,residuals] = regress(X_with_intercept, luminance(:));

% Now, residuals contain the part of Y that is not explained by X

% Assuming Z is another independent variable you're interested in
% Add a column of ones to Z to include an intercept in the model
Y_with_intercept = [rating(:)];

% Perform regression of residuals on Z
[beta2,beta2_int,residuals2] = regress(residuals, Y_with_intercept);
mdl = fitlm( residuals, Y_with_intercept, 'Intercept', false );
fprintf( '\n Residual beta p value: %0.4f', mdl.Coefficients.pValue )

% beta2 contains the regression coefficients of your model without the influence of X

pred_x = beta * luminance(:);
resid_x = pupil_size(:) - pred_x;

figure(3); clf;

plot( 1:numel(pred_x), pred_x, 'r', 'displayname', 'predicted pupil size' ); hold on;
plot( 1:numel(pupil_size), pupil_size, 'b', 'displayname', 'actual pupil size' );
legend;
xlabel( 'luminance' );
ylabel( 'pupil size' );

figure(4); clf;

% subplot( 1, 4, 1 );
% plot( 1:numel(pred_x), pred_x, 'r', 'displayname', 'predicted pupil size' ); hold on;
% plot( 1:numel(pupil_size), pupil_size, 'b', 'displayname', 'actual pupil size' );
% legend;
% xlabel( 'luminance' );
% ylabel( 'pupil size' );
% 
% subplot( 1, 4, 2 );
% scatter( luminance, pupil_size ); hold on;
% plot( linspace(-3, 3), beta * linspace(-3, 3) );
% xlabel( 'luminance' );
% ylabel( 'pupil size' );
% 
% subplot( 1, 4, 3 );
% scatter( rating, resid_x ); hold on;
% plot( linspace(-3, 3), beta2 * linspace(-3, 3) );
% xlabel( 'rating' );
% ylabel( 'residual pupil size' );
% 
% subplot( 1, 4, 4 );
% plot( rating, 'r', 'displayname', 'rating' ); hold on;
% plot( resid_x, 'b', 'displayname', 'residual pupil size' ); hold on;
% legend;
% xlabel( 'rating' );
% ylabel( 'residual pupil size' );


plot( t(:, 1:end-10), rating, 'r', 'displayname', 'rating' ); hold on;
plot( t(:, 1:end-10), resid_x, 'b', 'displayname', 'residual pupil size' ); hold on;
legend;
xlabel( 'rating' );
ylabel( 'residual pupil size' );


%%

function t = compute_shot_traces(samps, vid_infos, keep_vars)

pupil_size = [];
video_t = [];
affil_aggr_ratings = [];
affil_aggr_ratings_err = [];

rest_vars = table();

for i = 1:size(samps, 1)
  edf_info = samps.edf_info{i};
  
  num_samps = numel( edf_info.video_time );
  
  pupil_size = [ pupil_size; edf_info.pupil_size(1:num_samps) ];
  video_t = [ video_t; edf_info.video_time(1:num_samps) ];

  if ( 1 )
    match_start = ...
      samps.search_start_stop(i, 1) >= vid_infos.start & ...
      samps.search_start_stop(i, 1) <= vid_infos.stop;
  else
    match_start = vid_infos.start == samps.search_start_stop(i, 1);
  end
  
  match_code = vid_infos.code == samps.Code{i};
  match_rating = match_code & match_start;
  
  if ( nnz(match_rating) == 0 )
%     assert( false );
    mean_rating = nan;
    err_rating = nan;
  else
    mean_rating = nanmean( vid_infos.affil_aggr_slider_value(match_rating) );
    err_rating = plotlabeled.nansem( vid_infos.affil_aggr_slider_value(match_rating) );
  end
  
  affil_aggr_ratings = [ affil_aggr_ratings; ones(num_samps, 1) * mean_rating ];
  affil_aggr_ratings_err = [ affil_aggr_ratings_err; ones(num_samps, 1) * err_rating ];
  rest_vars = [ rest_vars; repmat(samps(i, keep_vars), num_samps, 1) ];
end

cleaned_pup = pupil_size;
cleaned_pup(cleaned_pup == 0) = nan;
cleaned_pup = fillmissing( cleaned_pup, 'linear' );
smoothed_pup = smoothdata( cleaned_pup, 'gaussian', 5e3 );
% smoothed_pup = smoothdata( cleaned_pup, 'gaussian', 1e3 );

t = table();
t.cleaned_pup = cleaned_pup;
t.smoothed_pup = smoothed_pup;
t.pupil_size = pupil_size;
t.video_t = video_t;
t.affil_aggr_ratings = affil_aggr_ratings;
t.affil_aggr_ratings_err = affil_aggr_ratings_err;
t = [ t, rest_vars ];

end

function t = compute_traces(samps, vid_infos, keep_vars)

pupil_size = [];
position = [];
video_t = [];
affil_aggr_ratings = [];
affil_aggr_ratings_err = [];

rest_vars = table();

for i = 1:size(samps, 1)
  edf_info = samps.edf_info{i};
  
  num_samps = numel( edf_info.video_time );
  
  pupil_size = [ pupil_size; edf_info.pupil_size(1:num_samps) ];
  video_t = [ video_t; edf_info.video_time(1:num_samps) ];
  position = [ position; edf_info.position(1:num_samps, :) ];

  if ( 1 )
    match_start = ...
      samps.search_start_stop(i, 1) >= vid_infos.start & ...
      samps.search_start_stop(i, 1) <= vid_infos.stop;
  else
    match_start = vid_infos.start == samps.search_start_stop(i, 1);
  end
  
  match_code = vid_infos.code == samps.Code{i};
  match_rating = match_code & match_start;
  
  if ( nnz(match_rating) == 0 )
%     assert( false );
    mean_rating = nan;
    err_rating = nan;
  else
    mean_rating = nanmean( vid_infos.affil_aggr_slider_value(match_rating) );
    err_rating = plotlabeled.nansem( vid_infos.affil_aggr_slider_value(match_rating) );
  end
  
  affil_aggr_ratings = [ affil_aggr_ratings; ones(num_samps, 1) * mean_rating ];
  affil_aggr_ratings_err = [ affil_aggr_ratings_err; ones(num_samps, 1) * err_rating ];
  rest_vars = [ rest_vars; repmat(samps(i, keep_vars), num_samps, 1) ];
end

cleaned_pup = pupil_size;
cleaned_pup(cleaned_pup == 0) = nan;
cleaned_pup = fillmissing( cleaned_pup, 'linear' );
smoothed_pup = smoothdata( cleaned_pup, 'gaussian', 5e3 );
% smoothed_pup = smoothdata( cleaned_pup, 'gaussian', 1e3 );

t = table();
t.cleaned_pup = cleaned_pup;
t.smoothed_pup = smoothed_pup;
t.pupil_size = pupil_size;
t.position = position;
t.video_t = video_t;
t.affil_aggr_ratings = affil_aggr_ratings;
t.affil_aggr_ratings_err = affil_aggr_ratings_err;
t = [ t, rest_vars ];

end

function traces = sort_traces(traces)

cols_to_sort = {'position', 'affil_aggr_ratings', 'affil_aggr_ratings_err', 'cleaned_pup', 'pupil_size', 'smoothed_pup', 'video_t'};

for i = 1:height(traces)
    time_col = traces(i, 'video_t');
    time_mat = cell2mat(table2array(time_col));
    [~, sorted_indices] = sort(time_mat);

    for p = 1:size(cols_to_sort, 2)
        col_name = cell2mat(cols_to_sort(p));
        col = traces(i, col_name);
        col_mat = cell2mat(table2array(col));
        if ( isvector(col_mat) )
          col_sorted = col_mat(sorted_indices);
        else
          col_sorted = col_mat(:, sorted_indices);
        end
        traces(i, col_name) = array2table({col_sorted});
    end

    % traces.smoothed_pup{i} = traces.smoothed_pup{i} - traces.smoothed_pup{i}(1);

end
end

function [boundaries, match_shot] = find_shot_boundaries(shots, time, code)

match_shot = find( ismember(shots.identifier, code) );

get_t = @(t) min( abs(time - t) );

boundaries = nan( numel(match_shot), 2 );
for i = 1:numel(match_shot)
  % for each shot ...
  t0 = shots.start(match_shot(i));
  t1 = shots.stop(match_shot(i));

  [~, ind0] = get_t( t0 );
  [~, ind1] = get_t( t1 );

  boundaries(i, :) = [ind0, ind1];
end

end