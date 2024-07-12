
raw = readtable( '/Users/Nick/Downloads/jamie-efe-pose-track/Fellini Videos 2024.csv' );
dates = cellfun( @datestr, strrep(raw{1, :}, 'Pilot Day ', ''), 'un', 0 );
video_files = raw{2:end, :};

tot_files = {};
tot_dates = {};
tot_indices = [];
for i = 1:size(video_files, 2)
  first_empty = find( cellfun(@isempty, video_files(:, i)), 1 );
  if ( isempty(first_empty) )
    first_empty = size(video_files, 1); 
  else
    first_empty = first_empty - 1;
  end
  curr_files = video_files(1:first_empty, i);
  vid_index = 1:first_empty;
  
  tot_files = [ tot_files; curr_files ];
  tot_dates = [ tot_dates; repmat(dates(i), numel(curr_files), 1) ];
  tot_indices = [ tot_indices; vid_index(:) ];
end

tot_files = strrep( tot_files, 'O', '0' );
vid_ids = strrep( tot_files, '.MP4', '' );
vid_info_table = table( ...
  string(tot_files), string(tot_dates), string(vid_ids) ...
  , tot_indices, 'va', {'file', 'date', 'id', 'index'} );

%%

% file_name = '/Users/Nick/Downloads/shorter_test_label_left_front.mp4.avi';
file_name = '/Users/Nick/Downloads/jamie-efe-pose-track/GX010025.avi';

vr = VideoReader( file_name );
frame_rate = vr.FrameRate;
if ( frame_rate == 0 )
  frame_rate = 24;
end

num_frames = max( 1, floor(vr.Duration * frame_rate) );
t = (0:num_frames-1) / frame_rate;

% rand_data = smoothdata( rand(num_frames, 1), 'gaussian', frame_rate * 2 );
% plot( t, rand_data );

%%

h5disp( from_file )

%%

% from_file = '~/Downloads/shorter_test_label_left_front.h5';
from_file = '/Users/Nick/Downloads/jamie-efe-pose-track/labels.v001.000_GX010025.analysis.h5';
from_files = shared_utils.io.find( ...
  '/Users/Nick/Downloads/jamie-efe-pose-track/monkey_tracking_output', '.h5' );
file_names = shared_utils.io.filenames( from_files );

[poses, tracks, mu, raw_tracks, node_names] = deal( {} );
for i = 1:numel(from_files)
  [poses{i}, tracks{i}, mu{i}, raw_tracks{i}, node_names{i}] = extract_keypoints( from_files{i} );
end

vid_index = cate1( ...
  arrayfun(@(i) ones(rows(poses{i}), 1) * i, 1:numel(poses), 'un', 0) );
vid_names = cate1(...
  arrayfun(@(i) strings(rows(poses{i}), 1) + file_names{i}, 1:numel(poses), 'un', 0) );

poses = vertcat( poses{:} );
tracks = vertcat( tracks{:} );
mu = vertcat( mu{:} );
raw_tracks = vertcat( raw_tracks{:} );
node_names = node_names{1};
vid_ids = extractBetween( vid_names, 17, 17 + numel('GX020026') - 1 );
[~, match_info] = ismember( vid_ids, vid_info_table.id );
pose_info_table = vid_info_table(match_info, :);

pose_info_table.vid_index = vid_index;
pose_info_table.vid_name = vid_names;
pose_info_table.manipulation(:) = "baseline";
pose_info_table.manipulation(contains(pose_info_table.date, '25-Jun-2024')) = 'psilo';

%%

left_ear_ind = find( strcmp(node_names, 'left_ear') );
left_eye_ind = find( strcmp(node_names, 'left_eye') );

right_ear_ind = find( strcmp(node_names, 'right_ear') );
right_eye_ind = find( strcmp(node_names, 'right_eye') );

nose_ind = find( strcmp(node_names, 'nose') );

relevant_node_inds = [ left_ear_ind, left_eye_ind, right_ear_ind, right_eye_ind, nose_ind ];
relevant_node_inds = left_ear_ind;
% relevant_node_inds = 1:numel(node_names);

% sel_poses = poses(frame_index, :);
sel_poses = poses;
sel_poses = reshape( sel_poses, size(sel_poses, 1), size(sel_poses, 2)/2, [] );
sel_poses = sel_poses + mu;
sel_poses = sel_poses(:, relevant_node_inds, :);

non_nan = find( ~any(isnan(sel_poses), [2, 3]) );
% non_nan = 1:size(sel_poses, 1);
% fi = non_nan(256);
fi = non_nan(1600);
sel_poses = sel_poses(fi, :, :);

face_center = squeeze( nanmean(sel_poses, 2) );

figure(1); clf; hold on;

frame = read( vr, fi );
imshow( frame, 'parent', gca );

colors = spring( numel(relevant_node_inds) );

for i = 1:size(sel_poses, 2)
  h = plot( gca, sel_poses(:, i, 1), sel_poses(:, i, 2), 'ko' );
  set( h, 'color', colors(i, :, :) );
end

h = plot( gca, face_center(1), face_center(2), 'ko' );
set( h, 'color', 'w', 'linewidth', 2 );

%%

instance_scores = h5read( from_file, '/instance_scores' );

%%

% overall movement
% limb rotation / movement
% head rotation

%%  variance in key points (aka, "scale" or "size" of body)

do_draw = false;

sizes = @(x, s) arrayfun(@(y) size(x, y), s);

centered_poses = raw_tracks - nanmean( raw_tracks, 2 );
plt_node_names = head_node_names();
[~, select_nodes] = ismember( plt_node_names, node_names );

matches_video = find( pose_info_table.file == "GX010025.MP4" );
matches_video = 1:size(pose_info_table, 1);

pose_scale = nan( size(raw_tracks, 1), 1 );

for i = 1:numel(matches_video)
  fprintf( '\n %d of %d', i, numel(matches_video) );
  
  ind = matches_video(i);
  pose_this_frame = squeeze( centered_poses(ind, select_nodes, :) );
  pose_to_plot = squeeze( raw_tracks(ind, select_nodes, :) );
  
  if ( do_draw )
    figure(1); clf();
    frame = read( vr, i );
    imshow( frame, 'parent', gca );
    hold on;
    gscatter( pose_to_plot(:, 1), pose_to_plot(:, 2), plt_node_names(:) );
  end

  [coeff, score, latent] = pca( pose_this_frame );
  
  if ( isempty(coeff) )
    continue
  end
  
  center = nanmean( pose_to_plot, 1 );
  latent = 2 * sqrt( latent );
  
  if ( do_draw )
    rot = acosd( coeff(1, 1) );  
    ellipse( latent(1), latent(2), rot, center(1), center(2) );
  end
    
  pose_scale(ind) = max( latent );
end

%%  estimate movement speed

rel_inds = 16;  % left ear;
% rel_inds = 1:size(raw_tracks, 2);

est_speed_on = raw_tracks;
if ( 1 )
  ref_pose_scale = nanmedian( pose_scale );
  est_speed_on = rescale_poses( est_speed_on, ref_pose_scale );
end
est_speed_on = est_speed_on(:, rel_inds, :);

gi = findeachv( vid_index );
tot_est_speed = zeros( rows(raw_tracks), 1 );
for i = 1:numel(gi)
  est_speed = estimate_movement_speed( est_speed_on(gi{i}, :, :), round(frame_rate * 0.1) );
  tot_est_speed(gi{i}) = est_speed;
end

%%

figure(2); clf;
[I, C] = findeach( pose_info_table, {'vid_index', 'vid_name', 'manipulation'} );

psilo_ind = C.manipulation == 'psilo';
base_ind = C.manipulation == 'baseline';
med_speed = rowifun( @nanmedian, I, tot_est_speed );
p = ranksum( med_speed(psilo_ind), med_speed(base_ind) );

[I, id, C] = rowsets( 4, C, {}, {}, 'manipulation', 'vid_name', 'to_string', 1 );
plots.simplest_barsets( med_speed, I, id, C ...
  , 'error_func', @plotlabeled.nansem ...
  , 'add_points', 1 ...
);

%%

figure(1); clf;

[I, C] = findeach( pose_info_table, {'index', 'manipulation', 'date'} );
med_speed = rowifun( @nanmedian, I, tot_est_speed );

[I, id, C] = rowsets( 3, C, 'manipulation', 'index', {}, 'to_string', 1 );
ord = plots.orderby( C, arrayfun(@num2str, 1:17, 'un', 0) );
[I, id, C] = rowref_many( ord, I, id, C );

plots.simplest_barsets( med_speed, I, id, C ...
  , 'as_line_plot', 1 ...
  , 'error_func', @plotlabeled.nansem ...
);

%%

figure(1); clf;
axs = plots.panels( [2, 1] );

overlay_on_video( axs(1), axs(2), vr, frame_rate, poses, mu, tot_est_speed ...
  , 'win_size', max(1, round(frame_rate * 10)) ...
  , 'smooth_size', max(1, round(frame_rate * 0.5)) ...
  , 'frame_offset', max(0, round(frame_rate * 60)) ...
);

%%

figure(1); clf;

poses = raw_tracks(:, rel_inds, :);
for i = 1:size(poses, 2)
  scatter( poses(:, i, 1), poses(:, i, 2) );
end

%%

figure(1); clf; hold on;

poses = raw_tracks(:, rel_inds, :);

for i = 1:size(poses, 2)
  gi = findeachv( vid_index );
  colors = hsv( numel(gi) );
  colors(end, :) = [0, 0, 0];
  
  for k = 1:numel(gi)
    xx = poses(gi{k}, i, 1);
    yy = poses(gi{k}, i, 2);
    
    scatter( xx, yy, 4, colors(k, :, :) );
  end
end

xlim( [0, 1920] );
ylim( [0, 1080] );

%%

function overlay_on_video(ax_image, ax_data, vr, frame_rate, poses, mus, data, varargin)

defaults = struct();
defaults.win_size = [];
defaults.smooth_size = [];
defaults.frame_offset = 0;

params = shared_utils.general.parsestruct( defaults, varargin );
assert( ~isempty(params.smooth_size) ...
  , 'expected provided (non-empty) value for smooth size' );
assert( ~isempty(params.win_size) ...
  , 'expected provided (non-empty) value for window size' );

frame_offset = params.frame_offset;

num_frames = max( 1, floor(vr.Duration * frame_rate) ) - frame_offset;

raw_data = nan( params.win_size, 1 );
next_index = 1;

for i = (frame_offset+1):num_frames
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
  
  rel_inds = 16;  % left ear
  
  sel_poses = poses(i, :);
  sel_poses = reshape( sel_poses, size(sel_poses, 1), size(sel_poses, 2)/2, [] );
  sel_poses = sel_poses + mus(i, :, :);
  centroid = squeeze( nanmean(sel_poses, 2) );
  
%   gscatter( sel_poses(sel_frame, :, 1)', sel_poses(sel_frame, :, 2)', good_node_names(:) );
  scatter( ax_image, sel_poses(1, rel_inds, 1)', sel_poses(1, rel_inds, 2)' );
  scatter( ax_image, centroid(1), centroid(2), 4, 'r*' );
  %%%
  
  plot( ax_data, smoothed_data );
  
  drawnow;
%   pause( 1/frame_rate );
  
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

function pts = head_node_names()

pts = {'right_eye', 'left_eye', 'nose', 'mouth' ...
  , 'right_ear', 'left_ear', 'neck_front', 'neck_back', 'headpost'};

end

function p = rescale_poses(raw_tracks, s)

centered = raw_tracks - nanmean( raw_tracks, 2 );
centered = centered .* s;
p = centered + nanmean( raw_tracks, 2 );

end

function [poses, tracks, mu, raw_tracks, good_node_names] = extract_keypoints(from_file)

tracks = h5read( from_file, '/tracks' );
node_names = deblank( h5read(from_file, '/node_names') );
good_node_names = { ...
    'right_elbow', 'right_hand', 'right_foot' ...
  , 'left_shoulder', 'left_elbow', 'left_hand', 'left_hip', 'left_knee', 'left_foot' ...
  , 'torso' ...
  , 'right_eye', 'left_eye', 'nose', 'mouth', 'right_ear', 'left_ear', 'neck_front', 'neck_back', 'headpost' ...
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