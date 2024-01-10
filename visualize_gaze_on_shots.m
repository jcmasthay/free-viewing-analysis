samps_p = 'D:\data\changlab\jamie\free-viewing\edf_samples';
samps_files = shared_utils.io.findmat( samps_p );
samps = shared_utils.io.fload( samps_files{3} );

vid_p = 'D:\data\changlab\jamie\free-viewing\videos';

%%

% function gen_shot(vid_p, shot)

pos_hist = zeros( 100, 2 );
pup_hist = zeros( size(pos_hist, 1), 1 );
frame_hist = zeros( size(pup_hist) );
is_fix = zeros( size(pup_hist) );
pos_ind = 1;

figure(1); clf; 
ax_im = subplot( 3, 1, 1 );
ax_pos = subplot( 3, 1, 2 );
ax_pup = subplot( 3, 1, 3 );

for idx = 1:size(samps, 1)
  
shot = samps(idx, :);

if ( strcmp(shot.block_type, 'C') )
  use_vid_p = fullfile( vid_p, 'scrambled' );
else
  use_vid_p = vid_p;
end

edf_info = shot.edf_info{1};
vr = VideoReader( fullfile(use_vid_p, shot.video_filename{1}) );

[I, frames] = findeachv( edf_info.video_frame );
for i = 1:numel(I)
  frame = read( vr, frames(i)+1 );
  bg_frame = zeros( [900, 1600, 3] );
  offi = floor( (size(bg_frame, 1:2 ) - size(frame, 1:2)) * 0.5 ) + 1;
  endi = offi + size( frame, 1:2 ) - 1;
  bg_frame(offi(1):endi(1), offi(2):endi(2), :) = double( frame ) / 255;
  
  pos = mean( edf_info.position(I{i}, :), 1 );
  pup = mean( edf_info.pupil_size(I{i}) );
  
  if ( pos_ind <= size(pos_hist, 1) )
    pos_hist(pos_ind, :) = pos(:)';
    pup_hist(pos_ind) = pup;
    frame_hist(pos_ind) = frames(i);
    pos_ind = pos_ind + 1;
  else
    pos_hist(1:end-1, :) = pos_hist(2:end, :);
    pos_hist(end, :) = pos;
    pup_hist(1:end-1) = pup_hist(2:end);
    pup_hist(end) = pup;
    frame_hist(1:end-1) = frame_hist(2:end);
    frame_hist(end) = frames(i);
  end
  
  is_fix(:) = 0;
  for j = 1:numel(edf_info.fixations)
    fix_frames = unique( edf_info.fixations(j).video_frame );
    is_fix(ismember(frame_hist, fix_frames)) = 1;
  end
  
  cla( ax_im ); cla( ax_pos ); cla( ax_pup ); hold( ax_im, 'on' );
  hold( ax_pos, 'on' );
  imshow( bg_frame, 'parent', ax_im );
  
  s = 50;
  rp = [ pos(1:2)-s*0.5, s, s ];
  if ( ~any(isnan(rp)) && all(pos < [1600, 900] & pos >= 0) )
    r = rectangle( ax_im, 'Position', rp );
    set( r, 'facecolor', 'r' );
  end
  
  tx = find( is_fix );
  
  plot( ax_pos, pos_hist );
  plot( ax_pup, pup_hist );
  scatter( ax_pos, tx, is_fix(is_fix > 0) * max(get(ax_pos, 'ylim')), 1 );
  title( ax_pup, 'pupil' ); 
  title( ax_pos, 'position' );
  drawnow;
%   pause( 1/60 );
end

end

% end