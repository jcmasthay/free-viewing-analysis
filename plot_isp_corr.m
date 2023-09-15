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

%%  inter scan path correlations

smooth_win_size = 100;  % ms, or [] for no smoothing
pos_samples = cellfun( @(x) x.position, edf_samples.edf_info, 'un', 0 );

if ( ~isempty(smooth_win_size) )
  do_smooth = @(x) smoothdata( x, 1, 'gaussian', smooth_win_size );
  pos_samples = cellfun( do_smooth, pos_samples, 'un', 0 );
end

corr_each = {'clip_id', 'start', 'block_type', 'affiliativeness', 'interactive_agency'};
mask = edf_samples.block_type == 'A';
mask(:) = true;
[I, isp_corr] = findeach( edf_samples, corr_each, mask );
isp_corr.num_seen = cellfun( @numel, I );

corr_coeffs = nan( numel(I), 2 );
for i = 1:numel(I)
  ind = I{i};
  if ( numel(ind) == 1 ), continue; end
  
  coeff_samples = [];  
  
  for j = 1:numel(ind)
    for k = 1:numel(ind)
      if ( j == k ), continue; end

      set0 = pos_samples{I{i}(j)};
      set1 = pos_samples{I{i}(k)};      
      n_corr = min( size(set0, 1), size(set1, 1) );      
      
%       num_sec_offset = 1;
%       offset = num_sec_offset * 1e3;
%       n_corr = n_corr - offset - 1;
%       
%       subset = offset:offset+n_corr;
%       if ( isempty(subset) )
%         continue
%       end
      subset = 1:min(size(set0, 1), size(set1, 1));
      
      rx = corr( set0(subset, 1), set1(subset, 1), 'rows', 'complete' );
      ry = corr( set0(subset, 2), set1(subset, 2), 'rows', 'complete' );
      
      coeff_samples(end+1, :) = abs( [rx, ry] );
    end
  end
  
  if ( isempty(coeff_samples) ), continue; end
  corr_coeffs(i, :) = nanmean( coeff_samples, 1 );
end

isp_corr.coeffs = corr_coeffs;

%%

if ( 0 )
  isp_corr = first_sec_only;
else
  isp_corr = all_secs;
end

plt_vec = nanmean( isp_corr.coeffs, 2 );

mask = isp_corr.block_type == 'A';
% mask = mask & isp_corr.affiliativeness ~= 'neutral';
% mask = mask & isp_corr.affiliativeness ~= 'neutral';
mask = mask & isp_corr.num_seen == 2;

[I, id, C] = rowsets( 4, isp_corr ...
  , {'block_type'} ...
  , {'interactive_agency'} ...
  , {'affiliativeness'} ...
  , {} ...
  , 'mask', mask, 'to_string', true );
C = strrep( C, '_', ' ');

figure(2); clf;
[axs, hs, xs] = plots.simplest_barsets( plt_vec, I, id, C ...
  , 'summary_func', @nanmedian ...
  , 'error_func', @plotlabeled.nansem ...
  , 'add_points', false ...
);

set( axs, 'xticklabelrotation', 10 );

%%

plt_vec = nanmean( isp_corr.coeffs, 2 );

mask = ~isnan( plt_vec );
mask = mask & isp_corr.affiliativeness ~= 'neutral';

[I, isp_stat_tbl] = findeach( isp_corr, {'block_type', 'interactive_agency'}, mask );
isp_stat_tbl.p = nan( rows(isp_stat_tbl), 1 );

for i = 1:numel(I)
  ind_a = intersect( find(isp_corr.affiliativeness == 'affiliative'), I{i} );
  ind_b = intersect( find(isp_corr.affiliativeness == 'aggressive'), I{i} );
  if ( ~isempty(ind_a) && ~isempty(ind_b) )
    isp_stat_tbl.p(i) = ranksum( plt_vec(ind_a), plt_vec(ind_b) );
  end
end

