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

%%  plot mean proportions of looking to rois, for clip types and block types

mask = look_props.affiliativeness ~= 'neutral';
mask = mask & look_props.interactive_agency ~= 'monkey_other_animal';

[I, id, C] = rowsets( 4, look_props ...
  , {} ...
  , {'affiliativeness', 'interactive_agency'} ...
  , 'block_type' ...
  , 'timestamp' ...
  , 'mask', mask, 'to_string', true );
C = strrep( C, '_', ' ');

figure(1); clf;
[axs, hs, xs] = plots.simplest_barsets( look_props.look_props, I, id, C ...
  , 'summary_func', @nanmean ...
  , 'error_func', @plotlabeled.nansem ...
  , 'add_points', false ...
);
ylim( axs, [0, 1] );
set( axs, 'xticklabelrotation', 10 );

%%  stat for above

group_names = {'affiliativeness', 'interactive_agency', 'block_type'};

props = look_props.look_props(mask);
prop_tbl = look_props(mask, :);
group = cellfun( @(x) ungroupi(findeach(prop_tbl, x)), group_names, 'un', 0 );

group_labels = fcat.from( cellstr(prop_tbl{:, group_names}), group_names );
anova_outs = dsp3.anovan( prop_tbl.look_props, group_labels, {}, group_names );

anova_outs.anova_tables{1}
anova_outs.comparison_tables{1}

