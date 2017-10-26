function [img, thresholded] = imagecleanup(filename, thresholdvalue)

if(nargin == 1)
    thresholdvalue = 150;
end

img = double(imread(filename));

bw = img(:, :, 1) + img(:, :, 2) + img(:, :, 3);

bwclosed = imclose(bw, strel('disk', 30));
img = bwclosed - bw;

bw2 = imdilate(img, strel('disk', 2));
img = bw2;

thresholded = bw2 < thresholdvalue;

imwrite(thresholded, sprintf('%s%s', filename, 'test.png'));

figure(1);
imagesc(thresholded); colormap(gray);

end