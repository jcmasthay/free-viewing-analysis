vr = VideoReader( '/Users/Nick/Downloads/shorter_test_label_left_front.mp4.avi' );
num_frames = max( 1, floor(vr.Duration * vr.FrameRate) );
t = (0:num_frames-1) / vr.FrameRate;

rand_data = smoothdata( rand(num_frames, 1), 'gaussian', vr.FrameRate * 2 );
plot( t, rand_data );

%%

from_file = '~/Downloads/shorter_test_label_left_front.h5';
[poses, tracks, mu, raw_tracks] = extract_keypoints( from_file );

%%

% overall movement
% limb rotation / movement
% head rotation

%%

est_speed = estimate_movement_speed( raw_tracks, round(vr.FrameRate * 0.1) );

%%

figure(1); clf; 
axs = plots.panels( [2, 1] );

overlay_on_video( axs(1), axs(2), vr, poses, mu, est_speed ...
  , 'win_size', max(1, round(vr.FrameRate * 10)) ...
  , 'smooth_size', max(1, round(vr.FrameRate * 0.5)) ...
);

%%

function overlay_on_video(ax_image, ax_data, vr, poses, mus, data, varargin)

defaults = struct();
defaults.win_size = [];
defaults.smooth_size = [];

params = shared_utils.general.parsestruct( defaults, varargin );
assert( ~isempty(params.smooth_size) ...
  , 'expected provided (non-empty) value for smooth size' );
assert( ~isempty(params.win_size) ...
  , 'expected provided (non-empty) value for window size' );

num_frames = max( 1, floor(vr.Duration * vr.FrameRate) );

raw_data = nan( params.win_size, 1 );
next_index = 1;

for i = 1:num_frames
  f = read( vr, i );
  
  if ( next_index <= numel(raw_data) )
    raw_data(next_index) = data(i);
    next_index = next_index + 1;
  else
    raw_data(1:end-1) = raw_data(2:end);
    raw_data(end) = data(i);
  end
  
  smoothed_data = smoothdata( raw_data, 'gaussian', params.smooth_size );
  
  cla( ax_image );  
  cla( ax_data );
  hold( ax_image, 'on' );
  hold( ax_data, 'on' );
  
  imshow( f, 'parent', ax_image );
  %%%
  
  sel_poses = poses(i, :);
  sel_poses = reshape( sel_poses, size(sel_poses, 1), size(sel_poses, 2)/2, [] );
  sel_poses = sel_poses + mus(i, :, :);
  centroid = squeeze( nanmean(sel_poses, 2) );
  
%   gscatter( sel_poses(sel_frame, :, 1)', sel_poses(sel_frame, :, 2)', good_node_names(:) );
  scatter( ax_image, sel_poses(1, :, 1)', sel_poses(1, :, 2)' );
  scatter( ax_image, centroid(1), centroid(2), 4, 'r*' );
  %%%
  
  plot( ax_data, smoothed_data );
  
  drawnow;
%   pause( 1/vr.FrameRate );
  
  if ( ptb.util.is_esc_down )
    break
  end
end

end

function est_speed = estimate_movement_speed(tracks, smooth_amt)

centroid = squeeze( nanmean(tracks, 2) );
grad_x = gradient( centroid(:, 1) );
grad_y = gradient( centroid(:, 2) );
est_speed = vecnorm( [grad_x, grad_y], 2, 2 );
est_speed = smoothdata( est_speed, 'gaussian', smooth_amt );

end

function [poses, tracks, mu, raw_tracks] = extract_keypoints(from_file)

tracks = h5read( from_file, '/tracks' );
node_names = deblank( h5read(from_file, '/node_names') );
good_node_names = { ...
    'right_elbow', 'right_hand', 'right_foot' ...
  , 'left_shoulder', 'left_elbow', 'left_hand', 'left_hip', 'left_knee', 'left_foot' ...
  , 'torso' ...
};

[~, lb] = ismember( good_node_names, node_names );

nan_inds = any( isnan(tracks), 3 );
p_missing = sum( nan_inds, 1 ) / size( nan_inds, 1 );
miss_frac = compose( "%s: %0.3f%%\n", string(node_names), 100*p_missing(:) );
[~, ord] = sort( p_missing );
miss_frac = miss_frac(ord);
fprintf( "Missing: %s", miss_frac );

tracks = tracks(:, lb, :);
mu = nanmean( tracks, 2 );
dev = nanstd( tracks, [], 2 );
raw_tracks = tracks;
tracks = tracks - mu;

sizes = @(x, d) arrayfun(@(s) size(x, s), d);
poses = reshape( tracks, [], prod(sizes(tracks, 2:3)) );

end