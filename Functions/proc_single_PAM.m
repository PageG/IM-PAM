function [Lmean, Lsd] = proc_single_PAM(image, n_leaves)
%proc_single_PAM Calculate average Fv'/Fm' for each sample (image)
%   Calculate Fv'/Fm' for light adapted leaf from a single measurement.
%   Output is a vector of average values per leaf.
%   Depends upon the 'im_pam_tiff_fvfm' and 'seg_leaf' functions.
%   Assumes PAR = 134 umol m-2 s-1 during data collection. Change as
%   appropriate.


if n_leaves == 1
    [FvFm] = im_pam_tiff_fvfm(image,134);
    
    for i = 1:size(FvFm,3)
        Lmean(i) = nanmean(nanmean(FvFm));
        Lsd(i) = std(reshape(FvFm, 1, 640*480), 'omitna');
    end
    
else
    [FvFm] = im_pam_tiff_fvfm(image,134);
    FvFms = seg_leaf(FvFm, FvFm, n_leaves);
    
    for i = 1:size(FvFms,3)
        Lmean(i) = nanmean(nanmean(FvFms(:,:,i)));
        Lsd(i) = std(reshape(FvFms(:,:,i), 1, 640*480), 'omitna');
    end
end



end

