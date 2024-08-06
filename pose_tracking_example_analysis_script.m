%%  load in the pose tracking data and video from one recording session

root_p = '/Users/Nick/Downloads/jamie-efe-pose-track';

vid_info = generate_pose_tracking_video_info_table( ...
  fullfile(root_p, 'Fellini Videos 2024.csv') );

% find all the video file paths (avi_ps) and names (avi_fs) in `root_p`
avi_ps = shared_utils.io.find( root_p, '.avi' );
avi_fs = string( shared_utils.io.filenames(avi_ps, false) );

% choose a video to examine
vid_ind = 1;
vr = VideoReader( avi_ps{vid_ind} );

% construct the file path to the pose tracking output file corresponding
% to the selected video
vid_id = avi_fs(vid_ind);
analysis_fname = compose( "labels.v001.000_%s.analysis.h5", vid_id );
analysis_h5_p = fullfile( root_p, 'monkey_tracking_output', analysis_fname );

% read the poses
[key_points, node_names] = extract_keypoints( analysis_h5_p );
s = compose( "frames: %d; key points: %d; components per key point: %d (x, y)" ...
  , size(key_points) );
fprintf( '\n%s\n', s );

%%  visualize poses

figure(1); clf; ax = gca();

% frame_ind = 128;
frame_ind = 256;
frame = read( vr, frame_ind );
pose = squeeze( key_points(frame_ind, :, :) );

show_pose_on_frame( ax, frame, pose, node_names );
show_face_vectors( ax, pose, node_names );

%%

thetas = compute_face_angles( key_points, node_names );
centroids = compute_centroids( key_points );
speeds = compute_speed( centroids );

%%

intruder_events_s = [ 60.2, 30.3, 128.7 ];
speed_traces = align_to_events( speeds, 30, intruder_events_s, -1, 5 );
  
%%

function speed = compute_speed(centroids)

validateattributes( centroids, {'double'}, {'2d'}, mfilename, 'cents' );
assert( size(centroids, 2) == 2, 'Expected centroids to be a matrix with 2 columns.' );

gx = gradient( centroids(:, 1) );
gy = gradient( centroids(:, 2) );
speed = vecnorm( [gx, gy], 2, 2 );

end

function cents = compute_centroids(key_points)

cents = squeeze( nanmean(key_points, 2) );

end

function thetas = compute_face_angles(key_points, node_names)

thetas = nan( size(key_points, 1), 1 );
for i = 1:size(key_points, 1)
  v0 = compute_face_vectors( squeeze(key_points(i, :, :)), node_names );
  thetas(i) = acos( v0(1) );
end

end

function [v0, v1, p0] = compute_face_vectors(pose, node_names)

validateattributes( pose, {'double'}, {'2d'}, mfilename, 'pose' );
assert( size(pose, 1) == numel(node_names) ...
  , 'Pose does not correspond to node names' );

le_ind = find( strcmp(node_names, 'left_eye') );
re_ind = find( strcmp(node_names, 'right_eye') );
nose_ind = find( strcmp(node_names, 'nose') );

assert( numel(le_ind) == 1 && numel(re_ind) == 1 && numel(nose_ind) == 1 ...
  , 'Missing at least one of "left_eye", "right_eye", "nose" in node names' );

p0 = pose(le_ind, :);
v0 = pose(re_ind, :) - p0;
v0 = v0 ./ norm( v0 );
v1 = [ -v0(2), v0(1) ];

end

function show_face_vectors(ax, pose, node_names)

[v0, v1, p0] = compute_face_vectors( pose, node_names );

line_len = 100;
p1 = p0 + v0 * line_len;
pm0 = mean( [p0(:), p1(:)], 2 );

h = quiver( ax, p0(1), p0(2), v0(1)*line_len, v0(2)*line_len );
set( h, 'linewidth', 1, 'color', 'r', 'displayname', 'face axis1' );

h = quiver( ax, pm0(1), pm0(2), v1(1)*line_len, v1(2)*line_len );
set( h, 'linewidth', 1, 'color', 'g', 'displayname', 'face axis2' );

rot = rad2deg( acos(v0(1)) );
text( ax, p0(1), p0(2), compose("Angle = %0.3f deg", rot), 'color', 'w' );

end

function show_pose_on_frame(ax, frame, pose, node_names)

axis( ax ); 
imshow( frame, 'parent', ax );
hold( ax, 'on' );
gscatter( pose(:, 1), pose(:, 2), strrep(node_names(:), '_', ' ') );

end

function [traces, t] = align_to_events(x, fps, events, look_back_t, look_ahead_t)

% align_to_events -- Align data to events.
%
%   [traces, t] = align_to_events( x, fps, events, look_back_t, look_ahead_t )
%   extracts a time-series drawn from samples of `x` separately for each 
%   event in `events`. Each time-series is a vector looking back `look_back_t`
%   and ahead `look_ahead_t` seconds relative to each event.
%
%   `traces` is a matrix with one row per event and with columns as
%   time-points. `t` identifies columns (timepoints) of `traces`.
%
%   `x` should be an array whose first dimension corresponds to time.

ind = max( 0, floor(events * fps) );
i0 = ind + round( look_back_t * fps );
i1 = ind + round( look_ahead_t * fps );

num_samples = round( look_ahead_t * fps ) - round( look_back_t * fps ) + 1;
traces = nan( numel(events), num_samples );

clns = colons( ndims(x)-1 );

for i = 1:numel(events)
  traces(i, :) = x((i0(i):i1(i))+1, clns{:});
end

t = round( look_back_t * fps ) : round( look_ahead_t * fps );
t = t ./ fps;

end

function [tracks, good_node_names] = extract_keypoints(from_file)

tracks = h5read( from_file, '/tracks' );

node_names = deblank( h5read(from_file, '/node_names') );

% only keep this subset of keypoints
good_node_names = { ...
    'right_elbow', 'right_hand', 'right_foot' ...
  , 'left_shoulder', 'left_elbow', 'left_hand', 'left_hip', 'left_knee', 'left_foot' ...
  , 'torso' ...
  , 'right_eye', 'left_eye', 'nose', 'mouth', 'right_ear', 'left_ear' ...
  , 'neck_front', 'neck_back', 'headpost' ...
};

[~, matched_node_indices] = ismember( good_node_names, node_names );
tracks = tracks(:, matched_node_indices, :);

end