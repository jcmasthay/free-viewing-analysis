data_p = '/Volumes/external3/data/changlab/jamie/free-viewing';
vid_name = 'Monkey Thieves S2E6.avi';

bbox_p = fullfile( data_p, 'detections', sprintf('%s-bbox', vid_name) );
vid_p = fullfile( data_p, 'videos' );

mats = shared_utils.io.findmat( bbox_p );
mat_names = shared_utils.io.filenames( mats );
frame_nums = cellfun( @(x) str2double(x(numel('bbox_')+1:end)), mat_names );

desired_frame = 800;

desired_ind = frame_nums == desired_frame;
bbox_file = load( mats{desired_ind} );

accept_detect = @(x) x.conf >= 0.1;
% [x_min, y_min, width, height]
bbox_file.detections = bbox_file.detections(cellfun(accept_detect, bbox_file.detections));

vid_reader = VideoReader( fullfile(vid_p, vid_name) );
frame = read( vid_reader, desired_frame + 1 );

figure( 1 ); clf; imshow( frame );
hold on;

im_dims = [ size(frame, 2), size(frame, 1) ];
draw_detections( gca, bbox_file.detections, im_dims );

%%

%%

function draw_detections(ax, detections, im_dims)

for i = 1:numel(detections)
  detect_rect = detections{i}.bbox;
  p0 = detect_rect(1:2) .* im_dims;
  wh = detect_rect(3:4) .* im_dims;
  rectangle( 'position', [p0, wh], 'Parent', ax );
end

end