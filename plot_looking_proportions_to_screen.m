samp_files = shared_utils.io.findmat( ...
  fullfile(fv_data_directory, 'edf_samples') ...
);

samp_sesh = shared_utils.io.filenames( samp_files );
samp_sesh = datestr(cellfun( @(x) x(1:numel('30-Aug-2023')), samp_sesh, 'un', 0) );
keep_files = datetime(samp_sesh) >= datetime('12012023', 'InputFormat', 'MMddyyyy');
samp_files = samp_files(keep_files);

vid_infos = shared_utils.io.fload( ...
  fullfile(fv_data_directory, 'videos/vid_info.mat') );

%%

scr_dims = [ 1600, 900 ];
vid_stats = table();

for i = 1:numel(samp_files)
  %%
  fprintf( '\n %d of %d', i, numel(samp_files) );
  
  samps = shared_utils.io.fload( samp_files{i} );
  vid_fnames = string( samps.video_filename );
  [~, match_vid] = ismember( vid_fnames, vid_infos.name );
  assert( all(match_vid > 0) );
  
  vid_size = vid_infos.size(match_vid, :);
  vid_rect = nan( size(vid_size, 1), 4 );
  
  for j = 1:size(vid_size, 1)
    vid_off = floor( max(0, scr_dims - vid_size(j, :)) * 0.5 ) + 1;
    vid_rect(j, :) = [ ...
      vid_off(1), vid_off(2), vid_off + vid_size(j, :) - 1 ];
  end
  
  p_ib = zeros( size(samps, 1), 1 );
  p_miss = zeros( size(p_ib) );
  dur_ib = zeros( size(p_ib) );
  
  for j = 1:size(samps, 1)
    edf_info = samps.edf_info{j};
    tf = shared_utils.rect.inside( ...
        vid_rect(j, :) ...
      , edf_info.position(:, 1) ...
      , edf_info.position(:, 2) ...
    );
    
    p_ib(j) = pnz( tf );
    dur_ib(j) = sum( tf );
    p_miss(j) = pnz( any(isnan(edf_info.position), 2) );
  end
  
  vid_stat = table( p_ib, dur_ib, p_miss, 'va' ...
    , {'prop_in_screen', 'duration_in_screen', 'prop_missing_data'} );
  
  vid_stats = [ vid_stats; [vid_stat, samps] ];
  
end

vid_stats.session = datetime(string(datestr(vid_stats.timestamp, 'mmddyy')) ...
  , 'InputFormat', 'MMddyy');

%%

% [mu_I, mu_tbl] = findeach( vid_stats, {'session', 'block_type'} );
% mu_tbl.prop_in_screen = cellfun( @(x) nanmean(vid_stats.prop_in_screen(x)), mu_I );

mu_tbl = vid_stats;
plt_var = 'prop_in_screen';
% plt_var = 'prop_missing_data';

[I, id, C] = rowsets( 3, mu_tbl ...
  , {}, {'session'}, {'block_type'} ...
  , 'to_string', 1 ...
);

[~, ord] = sort( datetime(datestr(C(:, 2))) );
[I, id, C] = rowref_many( ord, I, id, C );

figure(1); clf;
axs = plots.simplest_barsets( mu_tbl.(plt_var), I, id, C ...
  , 'as_line_plot', true ...
  , 'error_func', @plotlabeled.nansem ...
);

set( axs, 'xticklabelrotation', 45 );
shared_utils.plot.match_ylims( axs );
arrayfun( @(x) hold( x, 'on'), axs, 'un', 0 );
shared_utils.plot.add_vertical_lines( axs, 7.5 );