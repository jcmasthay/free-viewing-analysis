function [ps, matched] = match_files(src_ps, varargin)

for i = 1:numel(varargin)
  varargin{i} = cellstr( varargin{i} );
  assert( numel(varargin{i}) == 1, 'Expected char or scalar folder path.' );
end

src_ps = cellstr( src_ps );
ps = cell( numel(src_ps), numel(varargin) + 1 );
ps(:, 1) = src_ps;

matched = false( size(ps) );
matched(:, 1) = true;

src_fs = shared_utils.io.filenames( src_ps, true );

for i = 1:numel(varargin)
  v_ps = shared_utils.io.findmat( varargin{i} );
  v_fs = shared_utils.io.filenames( v_ps, true );
  [~, loc] = ismember( v_fs, src_fs );
  matched(:, i+1) = loc > 0;
  ps(:, i+1) = v_ps(loc(loc > 0));
end

end