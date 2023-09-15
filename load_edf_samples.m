function edf_samples = load_edf_samples(edf_sample_ps)

edf_samples = table();
for i = 1:numel(edf_sample_ps)
  edf_sample_f = shared_utils.io.fload( edf_sample_ps{i} );
  edf_info = edf_sample_f.edf_info;
  pos_var = cate1( cellfun(@(x) nanstd(x.position, [], 1), edf_info, 'un', 0) );
  edf_sample_f.position_variance = pos_var;
  edf_samples = [ edf_samples; edf_sample_f ];
end

edf_samples = convert_char_vars_to_string( edf_samples );

end