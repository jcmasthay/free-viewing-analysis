function merged = do_merge_shots(shot_table, start_var, stop_var, I, dur_thresh)

merged = table;

for i = 1:numel(I)
  fprintf( '\n %d of %d', i, numel(I) );
  
  subset = shot_table(I{i}, :);
  [~, ord] = sort( subset.(start_var) );
  subset = subset(ord, :);
  I{i} = I{i}(ord);

  accept_beg = 1;
  accept_end = 1;
  dur = subset.(stop_var) - subset.(start_var);

  while ( accept_beg <= rows(subset) && accept_end <= rows(subset) )
    if ( sum(dur(accept_beg:accept_end)) < dur_thresh )
      accept_end = accept_end + 1;
    else        
      first_subset = subset(accept_beg, :);
      first_subset.(stop_var) = subset.(stop_var)(accept_end);
      first_subset.merged_from = {I{i}(accept_beg:accept_end)};
      merged = [ merged; first_subset ];
      accept_beg = accept_end + 1;
      accept_end = accept_beg;
    end
  end
end

end