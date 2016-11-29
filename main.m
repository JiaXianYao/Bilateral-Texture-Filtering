im = im2double(imread('../Images/paper/angel.png'));
k = 5;
iter = 5;
J = bilateralTextureFilter(im, k, iter); figure; imshow(J);