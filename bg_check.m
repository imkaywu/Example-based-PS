%background check
addpath('/Users/BlacKay/Documents/Projects/Data/bps_src_data/Einstein');
mask_ref = imread('mask_ref_1.BMP');
mask_target = imread('mask_target_1.BMP');

edge_ref = edge(mask_ref, 'canny');
[center, radius, ~] = imfindcircles(edge_ref, [40, 90]);
edge_target = edge(mask_target, 'canny');
[edge_y, edge_x] = find(edge_target);

num_img = 26;
for i = 0 : num_img
    if i < 10
        im = imread(['img_1_0', num2str(i), '.BMP']);
    else
        im = imread(['img_1_', num2str(i), '.BMP']);
    end
    imshow(im); hold on;
    viscircles(center, radius, 'EdgeColor', 'r');
    plot(edge_x, edge_y, 'w.');
end