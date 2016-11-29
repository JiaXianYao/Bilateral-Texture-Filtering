% Bilateral filtering for color images
% Note: This function is based on the implementation by
% Douglas R. Lanman, Brown University
function B = bfilterColor(A, w, sigma_s, sigma_r)

    % Convert input RGB image to L*a*b color space
    A = applycform(A, makecform('srgb2lab'));
    
    % Rescale range variance
    sigma_r = 100 * sigma_r;
    
    % Pre-compute the Gaussian spatial kernel
    Gs = fspecial('gaussian', 2 * w + 1, sigma_s);

    % Create the waitbar
    h = waitbar(0, 'Applying bilateral filter...');
    set(h, 'Name', 'Bilateral Filter Progress');
    
    % Initialize the image
    dim = size(A);
    B = zeros(dim);
    
    for i = w + 1 : dim(1) - w
        for j = w + 1 : dim(2) - w
            
            % Extract the local patch
            P = A(i - w : i + w, j - w : j + w, :);
            
            % Compute the Gaussian intensity kernel
            dL = P(:, :, 1) - A(i, j, 1);
            da = P(:, :, 2) - A(i, j, 2);
            db = P(:, :, 3) - A(i, j, 3);
            Gr = exp(-(dL.^2+da.^2+db.^2)/(2 * sigma_r^2));
            
            % Compute the bilateral filter response
            F = Gs .* Gr;
            normF = sum(F(:));
            B(i, j, 1) = sum(sum(P(:, :, 1) .* F)) / normF;
            B(i, j, 2) = sum(sum(P(:, :, 2) .* F)) / normF;
            B(i, j, 3) = sum(sum(P(:, :, 3) .* F)) / normF;
            
        end
        waitbar(i / dim(1));
    end
    
    % Convert the filtered image back to RGB color space
    B = applycform(B, makecform('lab2srgb'));
    
    % Close the waitbar
    close(h);

end

