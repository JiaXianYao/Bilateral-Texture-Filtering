im = im2double(imread('../Images/woman-gray.png'));
im_gray = rgb2gray(im);


w = 3 * sigma_s;

gauss = fspecial('gaussian', 2 * w + 1, sigma_s);
im_gauss = conv2(im_gray, gauss, 'same');

im_filtered = bfilter(im_gray, w, sigma_s, sigma_r);

figure; imshow(im_gray);
figure; imshow(im_gauss);
figure; imshow(im_filtered);

% imshow(alpha, 'Colormap', jet(255));

im = im2double(imread('../Images/woman-gray.png'));
I = rgb2gray(im);
k = 7;
iter = 1;

tic;
J1 = bilateralTextureFilter(im, k, iter); figure; imshow(J1);
J2 = bilateralTextureFilterAlternative(im, k, iter); figure; imshow(J2);
toc;
