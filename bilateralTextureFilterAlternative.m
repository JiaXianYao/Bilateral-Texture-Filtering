% Bilateral texture filtering for images
% I: single or color channel image
% k: patch size (odd valued)
% iter: number of iterations
function J = bilateralTextureFilterAlternative(I, k, iter)

    % Check k is odd valued
    assert(mod(k, 2) == 1);
    
    % Parameters
    dimX = size(I, 1); % dimension of I in x
    dimY = size(I, 2); % dimension of I in y
    c = size(I, 3); % number of color channels
    s = 2 * k - 1; % spatial kernel size (2 * k - 1)
    half_s = floor(s / 2); % half spatial kernel size
    sigma_s = k - 1; % spatial sigma (k - 1)
    sigma_r = 0.025 * sqrt(c); % range sigma (0.05 * sqrt(c))
    
    % Initialize variables
    f = fspecial('gaussian', s, sigma_s); % Gaussian spatial kernel
    J = zeros(size(I));
    
    % Run the filter a few times
    for m = 1 : iter
        
        % Compute the grayscale image
        I_gray = rgb2gray(I);
        
        % Compute the blurred image
        B = boxBlur(I, k);
        
        % Compute the mRTV
        mRTV = computeMRTV(I_gray, k);
        
        % Compute the guidance image
        G_prime = computeGuidance(B, mRTV, k);
        
        % Compute the joint bilateral filtering        
        parfor i = 1 : dimX
            for j = 1 : dimY
                
                minX = max(i-half_s, 1);
                maxX = min(i+half_s, dimX);
                minY = max(j-half_s, 1);
                maxY = min(j+half_s, dimY);
                
                G_prime_patch = G_prime(minX:maxX, minY:maxY);
                g = exp(-(G_prime_patch - G_prime(i, j)).^2 ...
                    /(2 * sigma_r^2));
                fg = g .* f((minX:maxX)-i+half_s+1,(minY:maxY)-j+half_s+1);
                fg = fg ./ sum(fg(:));
                
                for k = 1 : c
                    I_patch = I(minX:maxX, minY:maxY, k);
                    J(i, j, k) = sum(I_patch(:) .* fg(:));
                end
                
            end
        end
        
        % Reset I for the new iteration
        I = J;

    end
    
end