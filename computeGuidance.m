% Compute the final guidance image (G') for bilateral texture filtering
% B:    blurred image (single or three channels)
% mRTV: computed mRTV values
% k:    patch size (odd valued)
function G_prime = computeGuidance(B, mRTV, k)

    % Check B and mRTV have the same dimension
    % and k is odd valued
    assert(size(B, 1) == size(mRTV, 1) && ...
           size(B, 2) == size(mRTV, 2) && mod(k, 2) == 1);
    
    % Parameters
    half_k = floor(k / 2); % half patch size
    dimX = size(B, 1); % dimension of B in x
    dimY = size(B, 2); % dimension of B in y
    sigma_alpha = 5*k; % used in equation (6)
    
    % Compute the initial guidance image with patch shift
    G = zeros(size(B));
    mRTV_min = zeros(size(mRTV));
    
    parfor i = 1 : dimX
        for j = 1 : dimY
            
            minX = max(1, i-half_k);
            minY = max(1, j-half_k);
            maxX = min(i+half_k, dimX);
            maxY = min(j+half_k, dimY);
            
            mRTV_patch = mRTV(minX:maxX, minY:maxY);
            mRTV_min(i, j) = min(mRTV_patch(:));
            
            [row, col] = find(mRTV_patch == mRTV_min(i, j), 1);
            G(i, j, :) = B(minX+row-1, minY+col-1, :);
            
        end
    end
    
    % Compute the alpha map in equation (6)
    alpha = 2 * ((1 ./ (1 + exp(-sigma_alpha * (mRTV - mRTV_min)))) - 0.5);
    
    % Compute the final guidance image
    G_prime = zeros(size(G));
    
    for i = 1 : size(G, 3)
        G_prime(:, :, i) = alpha .* G(:, :, i) + ...
            (1 - alpha) .* B(:, :, i);
    end
    
end