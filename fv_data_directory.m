function p = fv_data_directory(set_to)

persistent data_dir;

if ( isempty(data_dir) )
  data_dir = 'C:\source\data\free_viewing';
end

if ( nargin > 0 )
  data_dir = set_to;
end

p = data_dir;

end