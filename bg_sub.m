addpath('/Users/BlacKay/Documents/Projects/Data/bps_src_data/Einstein');
addpath('/Users/BlacKay/Documents/Projects/Open Source/livewire');

dir = '/Users/BlacKay/Documents/Projects/Data/bps_src_data/dinasaur/';
%% reference object
nviews = 2;
nimages = 35;
center = [124, 314]; %[123.5, 294.0]
radius = 103.5; % 95.3
for k = 1 : nviews
    for i = 1 : nimages
        if i <= 10
            name = [dir, 'img_', num2str(k - 1), '_0', num2str(i - 1), '.BMP'];
        else
            name = [dir, 'img_', num2str(k - 1), '_', num2str(i - 1), '.BMP'];
        end
        
        img = imread(name);
        imshow(img); hold on;
        viscircles(center, radius);
        hold off;
    end
    
    mask = imcircle(size(img, 2), size(img, 1), center, radius, 1);
    imwrite(mask, 'dire, mask.bmp');
end

%% target object
% img = imread([dir, 'img_1_00.BMP']);
% mask = livewire(img);
% imshow(mask);
% mask = 255 * uint8(mask);
% imwrite(mask, [dir, '/mask_tar_1.bmp']);
