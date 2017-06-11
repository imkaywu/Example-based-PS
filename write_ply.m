function write_ply(fname, h_map, c_map)

res_img = size(h_map, 1) * size(h_map, 2);
[y, x] = find(h_map);
ind = find(h_map);
num_pts = numel(x);

P = zeros(num_pts, 3);
C = zeros(num_pts, 3);
P(:, 1) = x;
P(:, 2) = y;
P(:, 3) = h_map(ind);
C(:, 1) = c_map(ind);
C(:, 2) = c_map(res_img + ind);
C(:, 3) = c_map(2 * res_img + ind);

header = 'ply\n';
header = [header, 'format ascii 1.0\n'];
header = [header, 'element vertex ', num2str(num_pts), '\n'];
header = [header, 'property float32 x\n'];
header = [header, 'property float32 y\n'];
header = [header, 'property float32 z\n'];
header = [header, 'property uchar red\n'];
header = [header, 'property uchar green\n'];
header = [header, 'property uchar blue\n'];
header = [header, 'end_header\n'];

data = [P, C];

fid = fopen(fname, 'w');
fprintf(fid, header);
dlmwrite(fname, data, '-append', 'delimiter', '\t', 'precision', 3);
fclose(fid);