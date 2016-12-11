% MAIN illustrates the use of BILATERALTEXTUREFILTER.
%
%    This demo shows the typical use of the bilateral texture filter
%    implemented in BILATERALTEXTUREFILTER. The function takes as input a
%    grayscale/color image I, patch size K, and number of iterations ITER.
%
%    Example images included in the "images/" folder have the following
%    naming convention. Each image is named "name-k-iter.ext", where "k"
%    and "iter" are the best patch size and number of iterations selected
%    for the image.
%
% JiaXian Yao, UC Berkeley, November 2016.



% Load an example image
I = im2double(imread('images/tiger-7-5.jpg'));

% Set patch size and number of iterations (listed in the image name)
k = 7;
iter = 5;

% Apply the bilateral texture filter
J = bilateralTextureFilter(I, k, iter);

% Display both the input and output images
figure;
set(gcf, 'Name', 'Bilateral Texture Filtering Result');

subplot(1,2,1); imshow(I);
title('Input Image');

subplot(1,2,2); imshow(J);
title('Result of Bilateral Texture Filtering');
