% Compute the final guidance image (G') based on Eq. (5) and (6)
%
% B:    blurred image (single or three color channeled)
% mRTV: computed mRTV values
% k:    odd-valued patch size

function G_prime = computeGuidance(B, mRTV, k)
    
    % Parameters
    dimX = size(B, 1);     % dimension of B in x
    dimY = size(B, 2);     % dimension of B in y
    half_k = floor(k / 2); % half of patch size
    sigma_alpha = 5*k;     % used in Eq. (6)
    
    % Compute the initial guidance image with patch shift
    G = zeros(size(B));
    mRTV_min = zeros(size(mRTV));
    
    parfor i = 1 : dimX
        for j = 1 : dimY
            
            % Extract the local patch
            minX = max(1, i-half_k);
            minY = max(1, j-half_k);
            maxX = min(i+half_k, dimX);
            maxY = min(j+half_k, dimY);
            mRTV_patch = mRTV(minX:maxX, minY:maxY);
            
            % Compute the guidance image
            mRTV_min(i, j) = min(mRTV_patch(:));
            [row, col] = find(mRTV_patch == mRTV_min(i, j), 1);
            G(i, j, :) = B(minX+row-1, minY+col-1, :);
            
        end
    end
    
    % Compute the alpha map in Eq. (6)
    alpha = 2 * ((1 ./ (1 + exp(-sigma_alpha * (mRTV - mRTV_min)))) - 0.5);
    
    % Compute the final guidance image
    G_prime = zeros(size(G));
    
    for i = 1 : size(G, 3)
        G_prime(:, :, i) = alpha .* G(:, :, i) + ...
            (1 - alpha) .* B(:, :, i);
    end
    
end