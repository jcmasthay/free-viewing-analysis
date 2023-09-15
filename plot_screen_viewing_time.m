%%  load edf samples

root_data_p = fv_data_directory;
edf_sample_fs = shared_utils.io.findmat( fullfile(root_data_p, 'edf_samples') );
edf_samples = load_edf_samples( edf_sample_fs );

%%

new_dend_table = load( '/Volumes/external3/data/changlab/jamie/free-viewing/data/dendro_table.mat' );

for i = 1:rows(edf_samples)
  clip_id = edf_samples.clip_id(i);
  match_id = find( new_dend_table.t.clip_id == clip_id );
  edf_samples.interactive_agency(i) = new_dend_table.t.interactive_agency(match_id);
end

%%

screen_rect = [ 0, 0, 1600, 900 ];
edf_infos = edf_samples.edf_info;

t_inside = zeros( numel(edf_infos), 1 );
p_inside = zeros( size(t_inside) );

for i = 1:numel(edf_infos)
  fprintf( '\n %d of %d', i, numel(edf_infos) );
  pos = edf_infos{i}.position;
  test_rect = @(x) shared_utils.rect.inside(screen_rect, pos(x, 1), pos(x, 2));
  ib_screen = arrayfun( test_rect, 1:rows(pos) );
  t_inside(i) = sum( ib_screen );
  p_inside(i) = t_inside(i) / sum( ~any(isnan(pos), 2) );
end

%%

plt_vec = p_inside;

mask = true( numel(plt_vec), 1 );

% mask = isp_corr.block_type == 'A';
% % mask = mask & isp_corr.affiliativeness ~= 'neutral';
% % mask = mask & isp_corr.affiliativeness ~= 'neutral';
% mask = mask & isp_corr.num_seen == 2;

[I, id, C] = rowsets( 4, edf_samples ...
  , {'block_type'} ...
  , {'interactive_agency'} ...
  , {'affiliativeness'} ...
  , {} ...
  , 'mask', mask, 'to_string', true );

C = strrep( C, '_', ' ');
C = strrep( C, 'A', 'chronological order block' );
C = strrep( C, 'B', 'shuffled order block' );
C = strrep( C, 'C', 'scrambled block' );

figure(2); clf;
[axs, hs, xs] = plots.simplest_barsets( plt_vec, I, id, C ...
  , 'summary_func', @nanmean ...
  , 'error_func', @plotlabeled.nansem ...
  , 'add_points', false ...
);

ylabel( axs(1), 'Proportion of clip' );
title( axs(1), 'chronological order block | proportion of time looking at screen' );

set( axs, 'xticklabelrotation', 10 );
ylim( axs, [0, 1.2] );