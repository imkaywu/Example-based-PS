% empty or filled circle

function circle = imcircle(w, h, c, r, t)
    if nargin == 4
        t = 0;
    end

    circle = zeros(h, w, 'uint8');
    x = repmat(1 : w, h, 1);
    y = repmat((1 : h)', 1, w);
    dist = sqrt((x - c(1)).^2 + (y - c(2)).^2);
    
    switch t
        case 0
            % empty circle
            ind_in_circ = find(dist <= r & dist > r - 0.5);
            circle(ind_in_circ) = 255;
            ind_out_circ = find(dist > r & dist < r + 0.5);
            circle(ind_out_circ) = uint8(255.0 * (r + 0.5 - dist(ind_out_circ)) / 0.5);
        case 1
            % filled circle
            ind_in_circ = find(dist <= r);
            circle(ind_in_circ) = 255;
            ind_on_circ = find(dist > r & dist - 0.5 < r);
            circle(ind_on_circ) = uint8(255.0 * (r + 0.5 - dist(ind_on_circ)) / 0.5);
    end
end