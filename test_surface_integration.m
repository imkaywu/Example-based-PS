clear;
close all;
[x, y] = meshgrid(-2:0.1: 2);
z = x.^2 + y.^2;
surf(x, y, z); hold on;
[nx, ny, nz] = surfnorm(x, y, z);
quiver3(x, y, z, nx, ny, nz);
nx = nx ./ nz;
ny = ny ./ nz;
mask = ones(size(x));

% Horn's method
Z = integrate_horn2(-nx, -ny, mask, 500, 1);

% Method 2
% Z = frankotchellappa(-nx, -ny);

% Shapelet method
% [slant, tilt] = grad2slanttilt(-nx, -ny);
% Z = shapeletsurf(slant, tilt, 6, 1, 2, 'slanttilt');
figure;
surf(x, y, Z);