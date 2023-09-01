function clip_table = append_clip_meta_data_to_clip_table(target_clips, clip_table, dend_table)

targ_clip_id = target_clips(clip_table.index, 'clip_id');
[~, loc] = ismember( targ_clip_id.clip_id, dend_table.clip_id );
clip_table = [ clip_table, dend_table(loc, :) ];

end