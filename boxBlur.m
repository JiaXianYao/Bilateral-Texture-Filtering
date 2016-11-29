% Box filter the input single-channel image
function B = boxBlur(I, k)

    % Parameters
    dimX = size(I, 1); % dimension of I in x
    dimY = size(I, 2); % dimension of I in y
    half_k = floor(k / 2); % half patch size

    % Initialize the image
    B = zeros(dimX, dimY);
    
    % Blur the image
    parfor i = 1 : dimX
        for j = 1 : dimY
            
            minX = max(i-half_k, 1);
            maxX = min(i+half_k, dimX);
            minY = max(j-half_k, 1);
            maxY = min(j+half_k, dimY);
            
            I_patch = I(minX:maxX, minY:maxY);
            B(i, j) = sum(I_patch(:)) / numel(I_patch);
            
        end
    end

end

