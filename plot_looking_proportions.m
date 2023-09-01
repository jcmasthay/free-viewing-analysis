%%  load preprocessed looking proportions

root_data_p = fv_data_directory();
look_prop_fs = shared_utils.io.findmat( fullfile(root_data_p, 'looking_proportions') );
look_props = cate1( cellfun(@shared_utils.io.fload, look_prop_fs, 'un', 0) );
look_props = convert_char_vars_to_string( look_props );

%%  plot histograms of clip types

[I, id, C] = rowsets( 1, look_props, {'affiliativeness', 'interactive_agency'} ...
  , 'to_string', true );
[PI, PL] = plots.nest1( id, I, C );

axs = plots.panels( numel(PI) );
for i = 1:numel(axs)
  plots.hist( axs(i), look_props.look_props(PI{i}{1}), strrep(PL{i, 1}, '_', ' ') );
  xlim( axs(i), [0, 1] );
end

shared_utils.plot.match_ylims( axs );

%%  proportions of looking to animal rois for affil vs. aggressive clip types

mask = look_props.interactive_agency == 'monkey_monkey';
% mask = true( rows(look_props), 1 );
mask = mask & look_props.block_type == 'A';

affil_ind = mask & look_props.affiliativeness == 'affiliative';
agg_ind = mask & look_props.affiliativeness == 'aggressive';

[p, ~, stats] = ranksum( look_props.look_props(affil_ind), look_props.look_props(agg_ind) )

%%  