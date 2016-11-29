function [B] = bfilter(A, w, sigma_s, sigma_r)

    if size(A, 3) == 1
        B = bfilterGray(A, w, sigma_s, sigma_r);
    else
        B = bfilterColor(A, w, sigma_s, sigma_r);
    end

end

