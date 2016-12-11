% BILATERALFILTER Two dimensional bilateral filtering.
%
%    This function implements 2-D bilateral filtering using the method
%    outlined in:
%
%        C. Tomasi and R. Manduchi. Bilateral Filtering for Gray and Color
%        Images. In Proceedings of the IEEE International Conference on
%        Computer Vision, 1998. 
%
%    B = bfilter2(A, W, SIGMA_S, SIGMA_R) performs 2-D bilateral filtering
%    for color image A. A should be a double precision matrix with
%    normalized values in the closed interval [0,1]. The half-size of the
%    Gaussian bilateral filter window is defined by W. The spatial-domain
%    standard deviation is given by SIGMA_S and the range-domain standard
%    deviation is given by SIGMA_R.
%
% JiaXian Yao, UC Berkeley, November 2016 (modified based on Douglas R.
% Lanman, Brown University, September 2006).

function B = bilateralFilter(A, w, sigma_s, sigma_r)

    % Verify the input image A is valid
    if ~isfloat(A) || size(A, 3) ~= 3 || min(A(:)) < 0 || max(A(:)) > 1
        error(['Input image A must be a double precision matrix of ', ...
               'size MxNx3 on the closed interval [0,1].']);
    end
    
    % Verify bilateral filter window size
    if isempty(w) || numel(w) ~= 1 || w < 1
        w = 5;
    end
    w = ceil(w);
    
    % Verify spatial-domain standard deviations.
    if isempty(sigma_s) || numel(sigma_s) ~= 1 || sigma_s <= 0
       sigma_s = 3;
    end
    
    % Verify range-domain standard deviations.
    if isempty(sigma_r) || numel(sigma_r) ~= 1 || sigma_r <= 0
       sigma_s = 0.1;
    end

    % Convert input RGB image to L*a*b color space
    A = applycform(A, makecform('srgb2lab'));
    
    % Rescale range variance
    sigma_r = 100 * sigma_r;
    
    % Pre-compute the Gaussian spatial kernel
    Gs = fspecial('gaussian', 2 * w + 1, sigma_s);

    % Parameters
    dimX = size(A, 1); % dimension of A in x
    dimY = size(A, 2); % dimension of A in y
    c = size(A, 3); % number of color channels of A
    
    % Initialize the image
    B = zeros(size(A));
    
    parfor i = 1 : dimX
        for j = 1 : dimY
            
            % Extract the local patch
            minX = max(i-w, 1);
            maxX = min(i+w, dimX);
            minY = max(j-w, 1);
            maxY = min(j+w, dimY);
            P = A(minX:maxX, minY:maxY, :);
            
            % Compute the Gaussian range kernel
            dL = P(:, :, 1) - A(i, j, 1);
            da = P(:, :, 2) - A(i, j, 2);
            db = P(:, :, 3) - A(i, j, 3);
            Gr = exp(-(dL.^2+da.^2+db.^2)/(2 * sigma_r^2));
            
            % Compute the bilateral filter response
            F = Gr .* Gs((minX:maxX)-i+w+1,(minY:maxY)-j+w+1);
            normF = sum(F(:));
            for k = 1 : c
                B(i, j, k) = sum(sum(P(:, :, k) .* F)) / normF;
            end
            
        end
    end
    
    % Convert the filtered image back to RGB color space
    B = applycform(B, makecform('lab2srgb'));

end