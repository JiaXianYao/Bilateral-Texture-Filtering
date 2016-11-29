% BILATERALTEXTUREFILTER Bilateral texture filtering.
%
%    This function implements bilateral texture filtering using the
%    algorithm described in:
%
%        Cho, H., Lee, H., Kang, H., and Lee, S. 2014. Bilateral texture
%        filtering. ACM Trans. Graphics 33, 4, 128:1-128:8
%
%    J = bilateralTextureFilter(I, k, iter) performs bilateral texture
%    filtering for the grayscale or color image I. I should be a double
%    precision matrix with normalized values in the closed interval [0,1].
%    The patch size is defined by K and it must be odd valued. The number
%    of filtering iterations is defined by ITER. All notations used in the
%    function are consistent with those used in the paper.
%
% JiaXian Yao, UC Berkeley, November 2016.

function J = bilateralTextureFilter(I, k, iter)

    % Check image and patch size
    assert(isfloat(I) && min(I(:)) >= 0 && max(I(:)) <= 1 ...
        && mod(k, 2) == 1);
    
    % Image parameters
    dimX = size(I, 1); % dimension of I in x
    dimY = size(I, 2); % dimension of I in y
    c = size(I, 3);    % number of color channels of I
    
    % Joint bilateral filtering parameters
    s = 2 * k - 1;             % spatial kernel size
    half_s = floor(s / 2);     % half of spatial kernel size
    sigma_s = k - 1;           % spatial sigma
    sigma_r = 0.025 * sqrt(c); % range sigma (0.05 * sqrt(c))
    
    % Pre-compute the Gaussian spatial kernel
    f = fspecial('gaussian', s, sigma_s);
    
    % Initialize matrices
    B = zeros(size(I));     % blurred image
    mRTVs = zeros(size(I)); % computed mRTV based on Eq. (4)
    J = zeros(size(I));     % output filtered image
    
    % Run the filtering a few times
    for m = 1 : iter
        
        % Compute the blurred image and mRTV
        for i = 1 : c
            [B(:, :, i), mRTVs(:, :, i)] = ...
                computeBlurAndMRTV(I(:, :, i), k);
        end
        mRTV = sum(mRTVs, 3) / c;
        
        % Compute the guidance image
        G_prime = computeGuidance(B, mRTV, k);
        
        % Compute the joint bilateral filtering        
        parfor i = 1 : dimX
            for j = 1 : dimY

                minX = max(i-half_s, 1);
                maxX = min(i+half_s, dimX);
                minY = max(j-half_s, 1);
                maxY = min(j+half_s, dimY);
                
                for k = 1 : c
                    
                    % Extract the local patch
                    G_prime_patch = G_prime(minX:maxX, minY:maxY, k);
                    
                    % Compute Gaussian range kernel
                    g = exp(-(G_prime_patch - G_prime(i, j, k)).^2 ...
                            /(2 * sigma_r^2));
                    fg = g .* f((minX:maxX)-i+half_s+1, ...
                        (minY:maxY)-j+half_s+1);
                    
                    % Compute the bilateral texture filter response
                    I_patch = I(minX:maxX, minY:maxY, k);
                    J(i, j, k) = sum(I_patch(:) .* fg(:)) / sum(fg(:));
                    
                end
                
            end
        end
        
        % Reset I for the new iteration
        I = J;

    end
    
end