load( '~/Downloads/data.mat' );
%%

samps_p = '/Volumes/external3/data/changlab/jamie/free-viewing/edf_samples';
bbox_p = '/Volumes/external3/data/changlab/jamie/free-viewing/detections';
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
sel

keep_vars = { 'identifier', 'block_type', 'Code', 'timestamp' };

all_traces = table();
for i = 1:numel(sel)
  traces = compute_traces( samps(I{sel(i)}, :), vid_infos, keep_vars );
  rest_vars = setdiff( traces.Properties.VariableNames, keep_vars );
  rest_traces = cellfun( @(x) {traces.(x)(:)'}, rest_vars );
  kept_vars = unique( traces(:, keep_vars) );
  assert( rows(kept_vars) == 1 );
  for j = 1:numel(rest_vars)
    kept_vars.(rest_vars{j}) = rest_traces(j);
  end
  all_traces = [ all_traces; kept_vars ];
end

all_traces = sort_traces(all_traces);

%%

codes = unique( all_traces(:, 'Code') );
for i = 1:size(codes, 1)
  bbox_mat = shared_utils.io.fload(...
    fullfile(bbox_p, sprintf('%s.avi-bbox', codes.Code{i}), 'all_bboxes.mat') );
  
  for j = 1:numel(bbox_mat)
    frame_bboxes = bbox_mat{j};
    for k = 1:numel(frame_bboxes)
      normalized_bbox = frame_bboxes{k}.bbox;
      scaled_bbox = normalized_bbox .* 
      
      error( 'xx' );
    end
  end
end

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

        lum = lum_table.luminance{string(lum_table.code) == all_traces.Code(ri(j))};
        lum_t = (0:numel(lum)-1)/fps;

        rounded_video_t = floor( video_t * round_to ) / round_to;

        for k = 1:numel(rounded_video_t)
            [~, ind] = min( abs(rounded_video_t(k) - time_axis) );
            aligned_traces(j, ind) = pupil_size(k); 
            aligned_ratings(j, ind) = ratings(k);
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
    end

    smoothed_lums = smoothdata( aligned_lums, 2, 'gaussian', 5e3 );
    smoothed_traces = smoothdata( aligned_traces, 2, 'gaussian', 5e3 );
    
    for p = 1:numel( ri )
        all_traces.aligned_lum( ri(p) ) = {smoothed_lums( p, : )};
        all_traces.aligned_trace( ri(p) ) = {smoothed_traces( p, : )};
        all_traces.aligned_rating( ri(p) ) = {aligned_ratings( p, : )};
        all_traces.aligned_time( ri(p) ) = {time_axis};
    end

end

%%

save_path = strcat('/Users/efedogruoz/Desktop/pupsize_plots/data.mat');
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

function traces = sort_traces(traces)

cols_to_sort = {'affil_aggr_ratings', 'affil_aggr_ratings_err', 'cleaned_pup', 'pupil_size', 'smoothed_pup', 'video_t'};

for i = 1:height(traces)
    time_col = traces(i, 'video_t');
    time_mat = cell2mat(table2array(time_col));
    [~, sorted_indices] = sort(time_mat);

    for p = 1:size(cols_to_sort, 2)
        col_name = cell2mat(cols_to_sort(p));
        col = traces(i, col_name);
        col_mat = cell2mat(table2array(col));
        col_sorted = col_mat(sorted_indices);
        traces(i, col_name) = array2table({col_sorted});
    end

    % traces.smoothed_pup{i} = traces.smoothed_pup{i} - traces.smoothed_pup{i}(1);

end
end