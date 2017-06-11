% Implementation of example-based photometric stereo by Hertzman and Seitz
close all; clear; clc;
addpath('/Users/BlacKay/Downloads/Software/flann-1.8.4-src/src/matlab');
%% read in images

num_images = 8;
imgs_tgt = cell(num_images, 1);
imgs_ref = cell(num_images, 1);

for i = 1 : num_images
   imgs_tgt{i} = imread(['photodata/greenbottle/bottle', num2str(i), '.png']);
   imgs_ref{i} = imread(['photodata/greenbottle/sphere', num2str(i), '.png']);
end

mask_tgt = imread('photodata/greenbottle/bottlemask.png');
mask_ref = imread('photodata/greenbottle/spheremask.png');

%% normal estimation

% OV vector
num_pts_tgt = sum(sum(mask_tgt ~= 0));
num_pts_ref = sum(sum(mask_ref ~= 0));

r_tgt = zeros(num_images, num_pts_tgt);
g_tgt = zeros(num_images, num_pts_tgt);
b_tgt = zeros(num_images, num_pts_tgt);
r_ref = zeros(num_images, num_pts_ref);
g_ref = zeros(num_images, num_pts_ref);
b_ref = zeros(num_images, num_pts_ref);

ind_tgt = find(mask_tgt ~= 0);
ind_ref = find(mask_ref ~= 0);
res_tgt = size(imgs_tgt{1}, 1) * size(imgs_tgt{1}, 2);
res_ref = size(imgs_ref{1}, 1) * size(imgs_ref{1}, 2);

for i = 1 : num_images
   r_tgt(i, :) = imgs_tgt{i}(ind_tgt)';
   g_tgt(i, :) = imgs_tgt{i}(res_tgt + ind_tgt)';
   b_tgt(i, :) = imgs_tgt{i}(2 * res_tgt + ind_tgt)';
   r_ref(i, :) = imgs_ref{i}(ind_ref)';
   g_ref(i, :) = imgs_ref{i}(res_ref + ind_ref)';
   b_ref(i, :) = imgs_ref{i}(2 * res_ref + ind_ref)';
end

ov_tgt = [r_tgt; g_tgt; b_tgt];
ov_ref = [r_ref; g_ref; b_ref];

% correspondence
load('index.mat');

% build_params.algorithm = 'kdtree';
% build_params.trees = ceil(log2(3 * num_images));
% [index, params, ~] = flann_build_index(ov_ref, build_params);
% 
% k = 1;
% search_params.checks = 3 * num_images;
% [ind_corr, dis_corr] = flann_search(index, ov_tgt, k, search_params);

% normal assignment
[y_tgt, x_tgt] = find(mask_tgt);
[y_ref, x_ref] = find(mask_ref);
ind_tgt = find(mask_tgt);
ind_ref = find(mask_ref);

y_min = min(y_ref);
y_max = max(y_ref);
x_min = min(x_ref);
x_max = max(x_ref);
y_center = (y_min + y_max) / 2.0;
x_center = (x_min + x_max) / 2.0;
radius = max(sqrt((y_ref - y_center).^2 + (x_ref - x_center).^2));

n_ref_vec = zeros(num_pts_ref, 3);
n_ref_vec(:, 1) = (x_ref - repmat(x_center, num_pts_ref, 1));
n_ref_vec(:, 2) = (y_ref - repmat(y_center, num_pts_ref, 1));
n_ref_vec(:, 3) = sqrt(repmat(radius + 2, num_pts_ref, 1).^2 - n_ref_vec(:, 1).^2 - n_ref_vec(:, 2).^2);
n_ref_vec = n_ref_vec ./ radius;
n_tgt_vec = n_ref_vec(ind_corr, :);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot the corresponding points
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dd = 5;
% imshow(imgs_ref{1});
% hold on;
% for r = min(y_ref) : dd : max(y_ref)
%     for c = min(x_ref) : dd : max(x_ref)
%         if mask_ref(r, c)
%             ind = find(ind_ref == sub2ind(size(mask_ref), r, c));
%             quiver3(c, r, 0, n_ref_vec(ind, 1), n_ref_vec(ind, 2), n_ref_vec(ind, 3), 5);
%         end
%     end
% end

% figure;imshow(imgs_tgt{1});
% hold on;
% for r = min(y_tgt) : dd : max(y_tgt)
%     for c = min(x_tgt) : dd : max(x_tgt)
%         if mask_tgt(r, c)
%             ind = find(ind_tgt == sub2ind(size(mask_tgt), r, c));
%             quiver3(c, r, 0, n_tgt_vec(ind, 1), n_tgt_vec(ind, 2), n_tgt_vec(ind, 3), 5);
%         end
%     end
% end

% subplot(1, 2, 1);
% imshow(imgs_tgt{1});
% hold on;
% subplot(1, 2, 2);
% imshow(imgs_ref{1});
% hold on;
% 
% for i = 191600 : 191700
%     subplot(1, 2, 1);
%     plot(x_tgt(i), y_tgt(i), 'r*');
% %     quiver3(x_tgt(i), y_tgt(i), 0, n_tgt_vec(i, 1), n_tgt_vec(i, 2), n_tgt_vec(i, 3), 0.5);
%     
%     subplot(1, 2, 2);
%     plot(x_ref(ind_corr(i)), y_ref(ind_corr(i)), 'r*');
% %     quiver3(x_ref(ind_corr(i)), y_ref(ind_corr(i)), 0, n_ref_vec(ind_corr(i), 1), n_ref_vec(ind_corr(i), 2), n_ref_vec(ind_corr(i), 3), 0.5);
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% n_ref_map = zeros([size(mask_ref), 3]);
% n_ref_map(ind_ref) = n_ref_vec(:, 1);
% n_ref_map(res_ref + ind_ref) = n_ref_vec(:, 2);
% n_ref_map(2 * res_ref + ind_ref) = n_ref_vec(:, 3);

n_tgt_map = zeros([size(mask_tgt), 3]);
n_tgt_map(ind_tgt) = n_tgt_vec(:, 1);
n_tgt_map(res_tgt + ind_tgt) = n_tgt_vec(:, 2);
n_tgt_map(2 * res_tgt + ind_tgt) = n_tgt_vec(:, 3);

save('normal_mask.mat', 'n_tgt_map', 'mask_tgt');


%% surface integration

% h_map = compute_heightMap(n_ref_map, mask_ref);
% h_map = compute_height_map(n_ref_map, mask_ref);
% [h_ref, w_ref] = size(mask_ref);
% figure;
% plot3(repmat(1 : w_ref, h_ref, 1), repmat((1 : h_ref)', 1, w_ref), h_map);

h_map = compute_heightMap(n_tgt_map, mask_tgt);
% h_map = compute_height_map(n_tgt_map, mask_tgt);
% [h_tgt, w_tgt] = size(mask_tgt);
% figure;
% plot3(repmat(1 : w_tgt, h_tgt, 1), repmat((1 : h_tgt)', 1, w_tgt), h_map);

write_ply('greenbottle.ply', h_map, imgs_tgt{1});