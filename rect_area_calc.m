function tot_area = rect_area_calc(rs)

base_area = sum( cellfun(@rect_area, rs) );
isect_area = 0;

for i = 1:numel(rs)-1
  for j = i+1:numel(rs)
    
    r0 = rs{i};
    r1 = rs{j};
  
    isect = [ ...
        max(r0(1), r1(1)) ...
      , max(r0(2), r1(2)) ...
      , min(r0(3), r1(3)) ...
      , min(r0(4), r1(4)) ...
    ];
    
    isect_area = isect_area + rect_area( isect );
  end
end

tot_area = base_area - isect_area;

end