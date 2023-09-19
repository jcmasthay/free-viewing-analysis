root_p = fv_data_directory;
dst_p = fullfile( root_p, 'clip_tables' );

shared_utils.io.require_dir( dst_p );
src_mats = shared_utils.io.findmat( fullfile(root_p, 'data'), true );
dend_table = shared_utils.io.fload( fullfile(root_p, 'data', 'dendro_table.mat') );

src_fname = shared_utils.io.filenames( src_mats, true );
accept_src_mat = cellfun( @(x) numel(x) == numel('11-Sep-2023 11_44_52.mat'), src_fname );
src_mats = src_mats(accept_src_mat);

for si = 1:numel(src_mats)
  fprintf( '\n %d of %d', si, numel(src_mats) );
  
  src_data = load( src_mats{si} );
  
  try
    clip_table = append_clip_meta_data_to_clip_table( ...
        src_data.file.params.target_clips ...
      , src_data.file.clip_table ...
      , dend_table ...
    );    
  
    src_p = src_mats{si};
    save( fullfile(dst_p, shared_utils.io.filenames(src_p, true)), 'clip_table' );
  catch err
    warning( err.message );
  end
end