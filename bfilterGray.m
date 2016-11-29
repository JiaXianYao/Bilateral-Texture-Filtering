% Bilateral filtering for grayscale images
% Note: This function is based on the implementation by
% Douglas R. Lanman, Brown University
function B = bfilterGray(A, w, sigma_s, sigma_r)

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
            P = A(i - w : i + w, j - w : j + w);
            
            % Compute the Gaussian intensity kernel
            Gr = exp(-(P - A(i, j)).^2/(2 * sigma_r^2));
            
            % Compute the bilateral filter response
            F = Gs .* Gr;
            B(i, j) = sum(P(:) .* F(:)) / sum(F(:));
            
        end
        waitbar(i / dim(1));
    end

    % Close the waitbar
    close(h);
    
end