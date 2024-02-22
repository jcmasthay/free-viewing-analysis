samps_p = '/Volumes/external3/data/changlab/jamie/free-viewing/edf_samples';
samps_files = shared_utils.io.findmat( samps_p );

%%  load survey data

infos_p = '~/Downloads';
infos = ["vid_infos.mat", "vid_infos2.mat", "vid_infos3.mat", "vid_infos4.mat", "vid_infos5.mat", "vid_infos6.mat"];
infos = cellfun( @(x) shared_utils.io.fload(fullfile(infos_p, x)), infos, 'un', 0 );
% first run of survey didn't randomize the initial value and didn't store
% this
infos{1}.affil_aggr_slider_initial_value(:) = nan;
for i = 1:3
  % adjust slider values when range was [0, 100] instead of [-50, 50]
  infos{i}.affil_aggr_slider_value = infos{i}.affil_aggr_slider_value - 50;
end
vid_infos = cate1( infos );

%%

vid_infos = readtable( '/Users/Nick/Downloads/all_data.xlsx' );
vid_infos.code = string( deblank(vid_infos.code) );
vid_infos.prolific_pid = string( deblank(vid_infos.prolific_pid) );
% vid_infos.timestamp = datetime( vid_infos.timestamp );

%%

lum_table = shared_utils.io.fload( '/Volumes/external3/data/changlab/jamie/free-viewing/videos/lum.mat' );

%%  search monkey samples for one clip

clip_id = 'Mcq_Human_Feed';
for i = 1:numel(samps_files)
  samps = shared_utils.io.fload( samps_files{i} );
  if ( any(contains(samps.Code, clip_id)) && all(strcmp(samps.block_type, 'A')) )
    disp( 'found this clip: ' );
    disp( i );
  end
end

%%  load one file

samps = shared_utils.io.fload( samps_files{123} );

%%  load all files

samps = cell( numel(samps_files), 1 );
for i = 1:numel(samps_files)
  fprintf( '\n %d of %d', i, numel(samps_files) );
  samps{i} = shared_utils.io.fload( samps_files{i} );
end
samps = cate1( samps );

%%  find shot times from survey data

samps.search_start_stop = get_merged_shot_starts_stops( samps );

%%  compute traces (e.g. pupil size) over course of movie watching, along with survey data

[I, C] = findeach( samps, {'timestamp', 'Code', 'block_type'} );
sel = find( ...
  strcmp(C.Code, 'Baboon_AggressiveDominanceDisplay') & ...
  strcmp(C.block_type, 'A'));

keep_vars = { 'identifier', 'Code', 'timestamp' };

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

%%

figure(1); clf;
axs = plots.panels( [3, 1] );
smooth_ax_ind = 1;
rating_ax_ind = smooth_ax_ind + 1;
lum_ax_ind = rating_ax_ind + 1;

title_str = plots.cellstr_join( ...
  {unique(all_traces(:, {'Code'}))} );
title_str = strjoin( title_str, ' | ' );
title_str = plots.strip_underscore( title_str );

ylim( axs(rating_ax_ind), [-50, 50] );

for i = 1:size(all_traces, 1)
  plot( axs(smooth_ax_ind), all_traces.video_t{i}, all_traces.smoothed_pup{i} );
  hold( axs(smooth_ax_ind), 'on' );
  
  plot( axs(rating_ax_ind), all_traces.video_t{i}, all_traces.affil_aggr_ratings{i} );
  hold( axs(rating_ax_ind), 'on' );
  
  if ( i == 1 )
    plot( axs(rating_ax_ind), all_traces.video_t{i}, all_traces.affil_aggr_ratings{i} - all_traces.affil_aggr_ratings_err{i} );
    plot( axs(rating_ax_ind), all_traces.video_t{i}, all_traces.affil_aggr_ratings{i} + all_traces.affil_aggr_ratings_err{i} );
  end
  
  lum = lum_table.luminance{string(lum_table.code) == all_traces.Code(i)};
  lum = smoothdata( lum, 'gaussian', floor(fps * 0.5) );
  
  fps = 30;
  plot( axs(lum_ax_ind), (0:numel(lum)-1)/fps, lum(:)' );
  
  title( axs(smooth_ax_ind), title_str );
  title( axs(rating_ax_ind), title_str );
  title( axs(lum_ax_ind), title_str );
  
  ylabel( axs(smooth_ax_ind), 'Pupil size' );
  ylabel( axs(rating_ax_ind), 'Affil vs. aggr rating' );
  ylabel( axs(lum_ax_ind), 'Video luminance' );
  
  shared_utils.plot.match_xlims( axs );
end

%%

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
    mean_rating = mean( vid_infos.affil_aggr_slider_value(match_rating) );
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