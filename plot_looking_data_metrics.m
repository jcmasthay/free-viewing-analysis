%%  load edf samples

root_data_p = fv_data_directory;
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

%%  

var_subset = edf_samples.position_variance(:, 1);

edf_samples.session = string( datestr(edf_samples.timestamp, 'mmddyyyy') );
[I, var_diff] = findeach( edf_samples, {'start', 'session', 'interactive_agency', 'affiliativeness'} );
var_diff.var = nan( numel(I), 1 );

for i = 1:numel(I)
  a_block = intersect( find(edf_samples.block_type == 'A'), I{i} );
  b_block = intersect( find(edf_samples.block_type == 'B'), I{i} );
  
  if ( numel(a_block) == 1 && numel(b_block) == 1 )
    var_diff.var(i) = var_subset(a_block) - var_subset(b_block);
  end
end

mask = var_diff.affiliativeness ~= 'neutral';

[I, id, C] = rowsets( 3, var_diff ...
  , {} ...
  , {'interactive_agency'} ...
  , {'affiliativeness'} ...
  , 'mask', mask, 'to_string', true );
C = strrep( C, '_', ' ');

figure(1); clf;
[axs, hs, xs] = plots.simplest_barsets( var_diff.var, I, id, C ...
  , 'summary_func', @nanmean ...
  , 'error_func', @plotlabeled.nansem ...
);

title( axs(1), 'Difference in variance in position for clips seen in vs. out of order' );

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

% title( 'Mean variance in position for clip and block types' );
ylabel( axs(1), 'Variance' );