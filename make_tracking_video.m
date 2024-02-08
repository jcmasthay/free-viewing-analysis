function make_tracking_video(vr, points, vw, varargin)

defaults = struct();
defaults.dst_frame_size_fraction = 1;

params = shared_utils.general.parsestruct( defaults, varargin );
colors = uint8( hsv(size(points, 2)) .* 255 );

num_frames = size( points, 1 );
for i = 1:num_frames
  fprintf( '\n %d of %d', i, num_frames );
  
  frame = read( vr, i );
  
  if ( params.dst_frame_size_fraction ~= 1 )
    frame = imresize( frame, params.dst_frame_size_fraction );
  end
  
  p = squeeze( points(i, :, :) );
%   p = p ./ [vr.Width, vr.Height ];
%   p = p ./ [1920, 1080];
  
  for j = 1:size(p, 1)
    frame = fill_rect( frame, p(j, :), [0.0125, 0.0125], colors(j, :) );
  end
  
  writeVideo( vw, frame );
end

end

function im = fill_rect(im, p, s, color)

if ( any(isnan(p)) )
  return
end

i0 = 1 + floor( size(im, 1) * p(2) );
i1 = i0 + max( 1, floor(size(im, 1) * s(2)) );

j0 = 1 + floor( size(im, 2) * p(1) );
j1 = j0 + max( 1, floor(size(im, 2) * s(1)) );

i1 = min( i1, size(im, 1) );
j1 = min( j1, size(im, 2) );

for i = i0:i1
  for j = j0:j1
    im(i, j, :) = color;
  end
end

end