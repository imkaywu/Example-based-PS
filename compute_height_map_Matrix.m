function h_map = compute_height_map_Matrix(n_map, mask)
    [h, w] = size(mask);
    res_map = w * h;
    
    [y, x] = find(mask);
    num_pts = numel(x);
    ind_map = find(mask);
    
    cols = [min(x); max(x)];
    rows = [min(y); max(y)];
    
    count = 1;
    seq_map = zeros(size(mask));
    for c = cols(1) : cols(2)
        ind = (x == c);
        seq_map(y(ind), c) = (count : count + sum(ind) - 1)';
        count = count + sum(ind);
    end
    
    row_empty = [];
    M = sparse(2 * num_pts, num_pts);
    u = sparse(2 * num_pts, 1);
    
    % if cols < rows
    for c = cols(1) : cols(2)
%         disp(c);
%         if(c == 252)
%         end
        
        ind_row = find(mask(:, c)); % the starting and ending ind_row
        
        seq_prev_col = seq_map(ind_row, c - 1);
        seq_curr_col = seq_map(ind_row, c);
        seq_next_col = seq_map(ind_row, c + 1);
        has_left_neighbor = seq_prev_col & seq_curr_col;
        has_right_neighbor = seq_curr_col & seq_next_col;
        
        has_bottom_neighbor = zeros(size(ind_row));
        has_top_neighbor = zeros(size(ind_row));
        has_bottom_neighbor(1 : end - 1) = (seq_curr_col(2 : end) - seq_curr_col(1 : end - 1) == 1); % 1: has a bottom neighbor, 0: has no bottom neighbor
        has_top_neighbor(2 : end) = (seq_curr_col(2 : end) - seq_curr_col(1 : end - 1) == 1); % 1: has a top neighbor, 0: has no top neighbor
        
        % 1st col: index of current pixels, 2nd: index of the right neighbors, 3rd: index of the bottom neighbors
        seq_neighbors = zeros(numel(ind_row), 3);
        seq_neighbors(:, 1) = seq_map(ind_row, c);
        seq_neighbors(has_right_neighbor(1 : end), 2) = seq_map(ind_row(has_right_neighbor(1 : end)), c + 1); % has right neighbor
        seq_neighbors(~has_right_neighbor(1 : end) & has_left_neighbor(1 : end), 2) = seq_map(ind_row(~has_right_neighbor(1 : end) & has_left_neighbor(1 : end)), c - 1); % has left neighbor, no right neighbor
        
        seq_neighbors(has_bottom_neighbor == 1, 3) = seq_neighbors(has_bottom_neighbor == 1, 1) + 1; % has bottom neighbor
        seq_neighbors(has_bottom_neighbor == 0 & has_top_neighbor == 1, 3) = seq_neighbors(has_bottom_neighbor == 0 & has_top_neighbor == 1, 1) - 1; % has top neighbor, no bottom neighbor
        
        % 1st col: left-right neighbor, 2nd col: up-bottom neighbor, 1: forward neighbor, 0, no neighbor, -1: backward neighbor
        num_neighbors = zeros(numel(ind_row), 2);
        num_neighbors(has_right_neighbor, 1) = 1;
        num_neighbors(~has_right_neighbor & has_left_neighbor, 1) = -1;
        num_neighbors(has_bottom_neighbor == 1, 2) = 1;
        num_neighbors(has_bottom_neighbor == 0 & has_top_neighbor == 1, 2) = -1;
        
        % forward x axis
        seq_pt = seq_neighbors(num_neighbors(:, 1) == 1, 1);
        seq_neighbor_pt = seq_neighbors(num_neighbors(:, 1) == 1, 2);
        u(2 * seq_pt - 1) = n_map(ind_map(seq_pt)); % n.x
        M(sub2ind(size(M), 2 * seq_pt - 1, seq_pt)) = n_map(2 * res_map + ind_map(seq_pt)); % n.z
        M(sub2ind(size(M), 2 * seq_pt - 1, seq_neighbor_pt)) = -n_map(2 * res_map + ind_map(seq_pt)); % n.z
        
        % forward y axis
        seq_pt = seq_neighbors(num_neighbors(:, 2) == 1, 1);
        seq_neighbor_pt = seq_neighbors(num_neighbors(:, 2) == 1, 3);
        u(2 * seq_pt) = n_map(res_map + ind_map(seq_pt)); % n.y
        M(sub2ind(size(M), 2 * seq_pt, seq_pt)) = n_map(2 * res_map + ind_map(seq_pt)); % n.z
        M(sub2ind(size(M), 2 * seq_pt, seq_neighbor_pt)) = -n_map(2 * res_map + ind_map(seq_pt));
        
        % backward x axis
        seq_pt = seq_neighbors(num_neighbors(:, 1) == -1, 1);
        seq_neighbor_pt = seq_neighbors(num_neighbors(:, 1) == -1, 2);
        u(2 * seq_pt - 1) = n_map(ind_map(seq_pt));
        M(sub2ind(size(M), 2 * seq_pt - 1, seq_pt)) = -n_map(2 * res_map + ind_map(seq_pt));
        M(sub2ind(size(M), 2 * seq_pt - 1, seq_neighbor_pt)) = n_map(2 * res_map + ind_map(seq_pt));
        
        % backward y axis
        seq_pt = seq_neighbors(num_neighbors(:, 2) == -1, 1);
        seq_neighbor_pt = seq_neighbors(num_neighbors(:, 2) == -1, 3);
        u(2 * seq_pt) = n_map(res_map + ind_map(seq_pt));
        M(sub2ind(size(M), 2 * seq_pt, seq_pt)) = -n_map(2 * res_map + ind_map(seq_pt));
        M(sub2ind(size(M), 2 * seq_pt, seq_neighbor_pt)) = n_map(2 * res_map + ind_map(seq_pt));
        
        % empty rows
        seq_pt = seq_neighbors(num_neighbors(:, 1) == 0, 1);
        if(~isempty(seq_pt))
            row_empty = [row_empty; 2 * seq_pt - 1];
        end
        seq_pt = seq_neighbors(num_neighbors(:, 2) == 0, 1);
        if(~isempty(seq_pt))
            row_empty = [row_empty; 2 * seq_pt];
        end
    end
    
    % if rows < cols, TBD
    
    % remove all empty rows
    M(row_empty, :) = [];
    u(row_empty, :) = [];
    
    % least squares
    h = (M.' * M) \ (M.' * u);
    
    % From sparse back to full matrix
    h = full(h);

    % Outliers due to singularity
    outlier_ind = abs(zscore(h))>10;
    h_min = min(h(~outlier_ind));
    h_max = max(h(~outlier_ind));

    % reassemble h back to 2D map
    count = 1;
    h_map = double(mask);
    for c = cols(1) : cols(2)
        ind = (x == c);
        h_map(y(ind), c) = (h(count : count + sum(ind) - 1) - h_min) / (h_max - h_min) * 255;
        count = count + sum(ind);
    end

end
