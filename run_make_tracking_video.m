h5_p = '/Users/Nick/Downloads/test_project_left_front_training_250frames_testlabelall.000_GX010020.analysis.h5';
% vid_p = '/Users/Nick/Downloads/GX010020.MP4';
vid_p = '/Users/Nick/Downloads/CS7632_Lectures/4_Path_Planning/7 - 2021_Game_AI_-_Search_Algorithms_-_Breadth_First_Search.mp4';
dst_p = fullfile( fileparts(vid_p), 'dst_vid.mp4' );

tracks = h5read( h5_p, '/tracks' );
vr = VideoReader( vid_p );
vw = VideoWriter( dst_p, 'MPEG-4' );
open( vw );

for i = 1:size(tracks, 1)
  tracks(i, :, :) = squeeze(tracks(i, :, :)) ./ [1920, 1080];
end

% tracks = tracks(1:1024, :, :);

make_tracking_video( vr, tracks, vw );

close( vw );