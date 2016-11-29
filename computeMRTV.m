% Compute mRTV for bilateral texture filtering
% I: single channel image
% k: patch size (odd valued)
function mRTV = computeMRTV(I, k)

    % Check image I is single channeled and k is odd valued
    assert(size(I, 3) == 1 && mod(k, 2) == 1);

    % Parameters
    eps = 10e-9; % epsilon used in equation (4)
    half_k = floor(k / 2); % half patch size
    dimX = size(I, 1); % dimension of I in x
    dimY = size(I, 2); % dimension of I in y
    
    % Compute the gradient magnitude
    Ixy = imgradient(I, 'sobel');
    
    % Compute mRTV for each pixel based on equation (4)
    mRTV = zeros(size(I));
    
    parfor i = 1 : dimX
        for j = 1 : dimY
            
            minX = max(i-half_k, 1);
            maxX = min(i+half_k, dimX);
            minY = max(j-half_k, 1);
            maxY = min(j+half_k, dimY);
            
            I_patch = I(minX:maxX, minY:maxY);
            Ixy_patch = Ixy(minX:maxX, minY:maxY);
            mRTV(i, j) = (max(I_patch(:)) - min(I_patch(:))) * ...
                max(Ixy_patch(:)) / (sum(Ixy_patch(:)) + eps);

        end
    end
    
end