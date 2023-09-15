function r = bbox_to_pixel_rect(bbox, im_dims, screen_dims)

p0 = bbox(1:2) .* im_dims;
wh = bbox(3:4) .* im_dims;

adj_x = (screen_dims(1) - im_dims(1)) * 0.5;
adj_y = (screen_dims(2) - im_dims(2)) * 0.5;

p0 = p0 + [ adj_x, adj_y ];
r = [ p0(1), p0(2), p0(1) + wh(1), p0(2) + wh(2) ];

end