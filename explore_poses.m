%%

% h5_f = '/Users/Nick/Downloads/test_project_left_front_training_250frames_testlabelall.000_GX010020.analysis.h5';
h5_f = '~/Downloads/shorter_test_label_left_front.h5';

tracks = h5read( h5_f, '/tracks' );
node_names = deblank( h5read(h5_f, '/node_names') );
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
tracks = tracks - mu;

sizes = @(x, d) arrayfun(@(s) size(x, s), d);
poses = reshape( tracks, [], prod(sizes(tracks, 2:3)) );

%%

non_nan = find( ~any(isnan(poses), 2) );

umap = UMAP();
umap_res = fit_transform( umap, poses(non_nan, :) );
umap_colors = jet( size(poses, 1) );
umap_colors = umap_colors(non_nan, :);

%%

nan_inds = find( any(isnan(poses), 2) );
non_nan = find( ~any(isnan(poses), 2) );

fprintf( 'Kept %0.3f%% (%d of %d) of frames\n' ...
  , numel(non_nan) / size(poses, 1) * 100, numel(non_nan), size(poses, 1) );

%%

tsne_res = tsne( poses(non_nan, :) );

base_colors = jet( size(poses, 1) );
tsne_colors = base_colors(non_nan, :);

%%

gmm_model = fitgmdist( tsne_res, 11, 'Replicates', 25 );

%%  estimate best K for k means clustering using silhouette score

num_iters = 100;

ks = 2:11;
ss = nan( num_iters, numel(ks) );

parfor it = 1:num_iters
  
fprintf( '\n %d of %d', it, num_iters );

s = nan( numel(ks), 1 );

ind = 1;
for k = ks
  idx = kmeans( tsne_res, k );
  s(ind) = mean( silhouette(tsne_res, idx) );
  ind = ind + 1;
end

ss(it, :) = s;

end

%%  plot silhouette scores

figure(1); clf;

[~, best_ks] = max( ss, [], 2 );
best_ks = ks(best_ks);
hist( best_ks, ks );
title( 'Proportion of iterations with highest silhouette score' );
xlabel( 'K' );
ylabel( 'Num iterations' );

cs = histcounts( best_ks, [ks, max(ks)+1] );
[~, ind] = max( cs );
best_k = ks(ind);

%%

[pred_label, cent] = kmeans( tsne_res, best_k );
% [pred_label, cent] = kmeans( poses(non_nan, :), 11 );
colors = jet( numel(unique(pred_label)) );
cluster_colors = colors(pred_label, :);

%%

[~, pred_label] = max( posterior(gmm_model, tsne_res), [], 2 );
colors = jet( numel(unique(pred_label)) );
cluster_colors = colors(pred_label, :);

[~, ind] = min( abs(tsne_res - [-27.3728, 34.5619]), [], 1 );

%%

res = tsne_res;
colors = tsne_colors;
% colors = cluster_colors;

% res = umap_res;
% colors = umap_colors;

mask = pred_label == 1;
mask(:) = true;

figure(2); clf;
scatter( res(mask, 1), res(mask, 2), 8, colors(mask, :) );

lim0 = min(min(get(gca, 'ylim')), min(get(gca, 'xlim')));
lim1 = max(max(get(gca, 'ylim')), max(get(gca, 'xlim')));
xlim( [lim0, lim1] );
ylim( [lim0, lim1] );

if ( 0 )
  hold on;
  voronoi( cent(:, 1), cent(:, 2) );
  text( cent(:, 1), cent(:, 2), compose("cluster %d", unique(pred_label)) );
end

%%  overlay points on frame for desired cluster

vr = VideoReader( '/Users/Nick/Downloads/shorter_test_label_left_front.mp4.avi' );

desired_cluster = 3;
cluster_ind = non_nan(pred_label == desired_cluster);

sel_frame = randsample( numel(cluster_ind), 1 );

sel_poses = poses(cluster_ind, :);
sel_poses = reshape( sel_poses, size(sel_poses, 1), numel(good_node_names), [] );
sel_poses = sel_poses + mu(non_nan(cluster_ind), :, :);

figure(1); clf;
imshow( read(vr, cluster_ind(sel_frame)) );

hold on;
gscatter( sel_poses(sel_frame, :, 1)', sel_poses(sel_frame, :, 2)', good_node_names(:) );

title( compose("Frame %d from cluster %d", sel_frame, desired_cluster) );

%%

desired_cluster = 15;
sampled_frame = randsample( find(pred_label == desired_cluster), 1 );

%%

pred_label(non_nan == 1150)

%%

figure(1); clf;

x_inds = non_nan(find(pred_label == 1));
bar( x_inds, ones(1, numel(x_inds)) );
xlim( [0, numel(pred_label) ] );

%%

vr = VideoReader( '/Users/Nick/Downloads/shorter_test_label_left_front.mp4.avi' );

figure(1); clf;
ax = gca;

sel_poses = reshape( poses, size(poses, 1), numel(good_node_names), [] );

pose_ind = 1;
for i = 1:size(poses, 1)
  fprintf( '\n %d of %d', i, size(poses, 1) );
  
  frame = read( vr, i );
  
  sel_poses = poses(i, :);
  sel_poses = reshape( sel_poses, size(sel_poses, 1), numel(good_node_names), [] );
  sel_poses = sel_poses + mu(i, :, :);
  
  colors = hsv( size(sel_poses, 2) );
  
  if ( 1 )
    cla( ax );
    hold( ax, 'on' );
    imshow( frame, 'parent', ax );
    scatter( ax, sel_poses(1, :, 1)', sel_poses(1, :, 2)', 10, colors );
  end
  
  f = getframe( gcf );
  
  if ( i == 1 )
    dst_vid_p = fullfile( vr.Path, 'overlaid.mp4' );
    vw = VideoWriter( dst_vid_p, 'MPEG-4' );
    open( vw );
  end
  
  writeVideo( vw, f );
end

close( vw );
