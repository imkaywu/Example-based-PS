function n_map_target = estimate_normal(data)

    addpath('/Users/BlacKay/Downloads/Software/flann-1.8.4-src/src/matlab');

%% read images
    num_img = data.num_img;
    
    mask_target = imread(data.name_mask_target);
    mask_target(mask_target > 0) = 1;
    mask_ref = imread(data.name_mask_ref);
    [center, radius, ~] = imfindcircles(mask_ref, round([0.3, 0.5] * min(size(mask_ref, 1), size(mask_ref, 2))));
    mask_ref(mask_ref > 0) = 1;
    
    img_target = cell(num_img, 1);
    img_ref = cell(num_img, 1);
    
    for i = 1 : data.num_img
        img = imread(data.name_img_target{i});
        size_img = size(img);
        img = reshape(img, [], 3);
        img(mask_target == 0, :) = 0;
        img_target{i} = reshape(img, size_img);
        
        img = imread(data.name_img_ref{i});
        size_img = size(img);
        img = reshape(img, [], 3);
        img(mask_ref == 0, :) = 0;
        img_ref{i} = reshape(img, size_img);
    end
    
%% construct OV
    size_mask_target = sum(sum(mask_target ~= 0));
    size_mask_ref = sum(sum(mask_ref ~= 0));

    r_target = zeros(num_img, size_mask_target);
    g_target = zeros(num_img, size_mask_target);
    b_target = zeros(num_img, size_mask_target);
    r_ref = zeros(num_img, size_mask_ref);
    g_ref = zeros(num_img, size_mask_ref);
    b_ref = zeros(num_img, size_mask_ref);

    ind_target = find(mask_target ~= 0);
    ind_ref = find(mask_ref ~= 0);
    res_target = size(img_target{1}, 1) * size(img_target{1}, 2);
    res_ref = size(img_ref{1}, 1) * size(img_ref{1}, 2);

    for i = 1 : num_img
       r_target(i, :) = img_target{i}(ind_target)';
       g_target(i, :) = img_target{i}(res_target + ind_target)';
       b_target(i, :) = img_target{i}(2 * res_target + ind_target)';
       r_ref(i, :) = img_ref{i}(ind_ref)';
       g_ref(i, :) = img_ref{i}(res_ref + ind_ref)';
       b_ref(i, :) = img_ref{i}(2 * res_ref + ind_ref)';
    end

    ov_ref = [r_ref; g_ref; b_ref];
    ov_target = [r_target; g_target; b_target];
    
%% correspondence
    build_params.algorithm = 'kdtree';
    build_params.trees = ceil(log2(3 * num_img));
    [index, params, ~] = flann_build_index(ov_ref, build_params);

    k = 1;
    search_params.checks = 3 * num_img;
    [ind_corr, dis_corr] = flann_search(index, ov_target, k, search_params);
    
%% normal estimation
    [y_ref, x_ref] = find(mask_ref);
    n_ref_vec = zeros(size_mask_ref, 3);
    n_ref_vec(:, 1) = x_ref - repmat(center(1), size_mask_ref, 1);
    n_ref_vec(:, 2) = y_ref - repmat(center(2), size_mask_ref, 1);
    
    r = max(sqrt(n_ref_vec(:, 1).^2 + n_ref_vec(:, 2).^2));
    if(radius < r)
        radius = r + 2; % why 2?
    end
    
    n_ref_vec(:, 3) = sqrt(repmat(radius, size_mask_ref, 1).^2 - n_ref_vec(:, 1).^2 - n_ref_vec(:, 2).^2);
    n_ref_vec = n_ref_vec ./ radius;
    n_target_vec = n_ref_vec(ind_corr, :);
    
    % reshape normal
    n_map_ref = zeros([size(mask_ref), 3]);
    n_map_ref(ind_ref) = n_ref_vec(:, 1);
    n_map_ref(res_ref + ind_ref) = n_ref_vec(:, 2);
    n_map_ref(2 * res_ref + ind_ref) = n_ref_vec(:, 3);
    
    n_map_target = zeros([size(mask_target), 3]);
    n_map_target(ind_target) = n_target_vec(:, 1);
    n_map_target(res_target + ind_target) = n_target_vec(:, 2);
    n_map_target(2 * res_target + ind_target) = n_target_vec(:, 3);
    
%% show surface normal
%     show_surfNorm(n_map_ref, 5);
%     show_surfNorm(n_map_target, 5);
    
%% surface integration
    n_z = n_map_ref(:, :, 3);
    n_z(n_z == 0) = 1;
    n_x = -n_map_ref(:, :, 1) ./ n_z;
    n_y = -n_map_ref(:, :, 2) ./ n_z;
    h_ref = integrate_horn2(n_x, n_y, double(mask_ref), 10000, 1);
    write_ply('ref.ply', h_ref, img_ref{1});
    
    n_z = n_map_target(:, :, 3);
    n_z(n_z == 0) = 1;
    n_x = -n_map_target(:, :, 1) ./ n_z;
    n_y = -n_map_target(:, :, 2) ./ n_z;
    h_target = integrate_horn2(n_x, n_y, double(mask_target), 10000, 1);
    write_ply('target.ply', h_target, img_target{1});
end