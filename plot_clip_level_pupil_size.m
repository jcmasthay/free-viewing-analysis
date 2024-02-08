samps_p = '/Volumes/external3/data/changlab/jamie/free-viewing/edf_samples';

samps_files = shared_utils.io.findmat( samps_p );
samps = shared_utils.io.fload( samps_files{75} );

% clip_id = 'Mcq_Ingroup_SettlementBattle1';
% for i = 1:numel(samps_files)
%   samps = shared_utils.io.fload( samps_files{i} );
%   if ( any(contains(samps.Code, clip_id)) )
%     disp( 'found this clip: ' );
%     disp( i );
%   end
% end

%%

pupil_size = [];
video_t = [];

for i = 1:size(samps, 1)
  edf_info = samps.edf_info{i};
  
  num_samps = numel( edf_info.video_time );
  
  pupil_size = [ pupil_size; edf_info.pupil_size(1:num_samps) ];
  video_t = [ video_t; edf_info.video_time(1:num_samps) ];
end

cleaned_pup = pupil_size;
cleaned_pup(cleaned_pup == 0) = nan;
cleaned_pup = fillmissing( cleaned_pup, 'linear' );

smoothed_pup = smoothdata( cleaned_pup, 'gaussian', 5e3 );

title_str = strjoin( unique(samps.Code), ' | ' );
title_str = strrep( title_str, '_', ' ');

figure(1); 
clf;

axs = plots.panels( [3, 1] );
plot( axs(1), video_t, pupil_size );
plot( axs(2), video_t, cleaned_pup );
plot( axs(3), video_t, smoothed_pup );
arrayfun( @(x) title(x, title_str), axs );