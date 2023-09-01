%%  load edf samples

edf_sample_fs = shared_utils.io.findmat( fullfile(root_data_p, 'edf_samples') );
edf_samples = table();
for i = 1:numel(edf_sample_fs)
  edf_sample_f = shared_utils.io.fload( edf_sample_fs{i} );
  edf_info = edf_sample_f.edf_info;
  pos_var = cate1( cellfun(@(x) nanstd(x.position, [], 1), edf_info, 'un', 0) );
  edf_sample_f.position_variance = pos_var;
  edf_samples = [ edf_samples; edf_sample_f ];
end

edf_samples = convert_char_vars_to_string( edf_samples );

%%  plot histograms of variance in position for clip types

mask = edf_samples.block_type == 'A';
[I, id, C] = rowsets( 1, edf_samples, {'affiliativeness', 'interactive_agency'} ...
  , 'to_string', true, 'mask', mask );
[PI, PL] = plots.nest1( id, I, C );

plot_var = edf_samples.position_variance(:, 1);

axs = plots.panels( numel(PI) );
for i = 1:numel(axs)
  plots.hist( axs(i), plot_var(PI{i}{1}), strrep(PL{i, 1}, '_', ' ') );
end

shared_utils.plot.match_ylims( axs );

%%  plot means of variance in position for clip types and block types

mask = edf_samples.affiliativeness ~= 'neutral';

[I, id, C] = rowsets( 3, edf_samples ...
  , {'interactive_agency'} ...
  , {'affiliativeness'} ...
  , 'block_type' ...
  , 'mask', mask, 'to_string', true );
C = strrep( C, '_', ' ');

plot_var = edf_samples.position_variance(:, 1);

figure(1); clf;
[axs, hs, xs] = plots.simplest_barsets( plot_var, I, id, C ...
  , 'summary_func', @nanmean ...
  , 'error_func', @plotlabeled.nansem ...
);