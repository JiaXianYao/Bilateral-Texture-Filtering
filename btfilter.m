% Bilateral texture filtering for grayscale images
% I: input image
% k: patch size (odd valued)
function J = btfilter(I, k)

    % Parameters
    eps = 10e-9;
    hk = floor(k / 2); % half patch size
    dim = size(I); % image size
    sigmaA = 5*k;
    s = 2 * k - 1; % spatial kernel size
    hs = floor(s / 2); % half spatial kernel size
    sigmaS = k - 1;
    sigmaR = 0.05;

    % Uniform blurring of I with a box kernel
    box = ones(k, k) / k^2;
    B = conv2(I, box, 'same');
    
    % Compute mRTV for each pixel (with a waitbar)
    Ix = conv2(I, [-1, 0, 1], 'same');
    Iy = conv2(I, [-1; 0; 1], 'same');
    Ixy = Ix.*Ix + Iy.*Iy;
    Ixy = Ixy.^0.5;
    
    mRTV = zeros(size(I));
    h = waitbar(0, 'Computing mRTV...');
    set(h, 'Name', 'Progress');
    
    for i = 1 : dim(1)
        for j = 1 : dim(2)
            
            minX = max(1, i-hk);
            minY = max(1, j-hk);
            maxX = min(i+hk, dim(1));
            maxY = min(j+hk, dim(2));
            
            pI = I(minX:maxX, minY:maxY);
            pIxy = Ixy(minX:maxX, minY:maxY);
            mRTV(i, j) = (max(pI(:)) - min(pI(:))) * ...
                max(pIxy(:)) / (sum(pIxy(:)) + eps);
            
        end
        waitbar(i / dim(1));
    end
    
    close(h);
    
    % Compute the guidance image with patch shift (with a waitbar)
    G = zeros(size(I));
    mRTVMin = zeros(size(I));
    h = waitbar(0, 'Computing guidance image...');
    set(h, 'Name', 'Progress');
    
    for i = 1 : dim(1)
        for j = 1 : dim(2)
            
            minX = max(1, i-hk);
            minY = max(1, j-hk);
            maxX = min(i+hk, dim(1));
            maxY = min(j+hk, dim(2));
            
            p = mRTV(minX:maxX, minY:maxY);
            pB = B(minX:maxX, minY:maxY);
            mRTVMin(i, j) = min(p(:));
            G(i, j) = pB(p == mRTVMin(i, j));
            
        end
        waitbar(i / dim(1));
    end
    
    close(h);
    
    % Compute alpha and the new guidance image
    alphas = 2 * ((1 ./ (1 + exp(-sigmaA * (mRTV - mRTVMin)))) - 0.5);
    Gp = alphas .* G + (1 - alphas) .* B;
    
    % Appy joint bilateral filtering
    
    % Pre-compute the Gaussian spatial kernel
    f = fspecial('gaussian', s, sigmaS);

    % Create the waitbar
    h = waitbar(0, 'Applying joint bilateral filter...');
    set(h, 'Name', 'Joint Bilateral Filter Progress');
    
    % Initialize the image
    J = zeros(dim);
    
    for i = 1 : dim(1)
        for j = 1 : dim(2)
            
            % Extract the local patch
            minX = max(1, i-hs);
            minY = max(1, j-hs);
            maxX = min(i+hs, dim(1));
            maxY = min(j+hs, dim(2));
            p = Gp(minX:maxX, minY:maxY);
            Ip = I(minX:maxX, minY:maxY);
            
            % Compute the range kernel
            g = exp(-(p - Gp(i, j)).^2/(2 * sigmaR^2));
            
            % Compute the bilateral filter response
            fg = g .* f(hs+1-(i-minX):hs+1+(maxX-i), ...
                        hs+1-(j-minY):hs+1+(maxY-j));
            J(i, j) = sum(Ip(:) .* fg(:)) / sum(fg(:));
            
        end
        waitbar(i / dim(1));
    end

    % Close the waitbar
    close(h);
    
end

