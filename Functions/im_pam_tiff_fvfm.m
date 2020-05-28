function [FvFm] = im_pam_tiff_fvfm(file,PAR)
%im_pam_tiff_fvfm Import TIFF stack from Waltz Imaging PAM
%   Use this function to import and process and TIFF stack exported from
%   ImagingWin software (tested with ImagingWin v2.46i, with data collected using an M-series Imaging PAM
%   fitted with an IMAG-MIN/B measuring head, Walz GmbH, Effeltrich, Germany. Outputs are calculated from
%   formalae described in the Imaging PAM manual.
%   PAR is a single value.
%   This function only calculates Fv'Fm' for single measurements. If more
%   parameters are required, see 'im_pam_tiff.m'

fname = file;
info = imfinfo(fname);
num_images = numel(info);

flstack = zeros(480, 640, num_images);

for k = 1:num_images
    A = imread(fname, k, 'Info', info);
    A = double(A);
    A(A < 10) = NaN; % filter out noise
    flstack(:,:,k) = A;
    A(A == 0) = NaN;
    flstack(:,:,k) = flstack(:,:,k)*10^-3;
    clear A
end

% Fo, Fm, NIR, Red, F1, Fm'1, F2, Fm'2, F3, Fm'3 etc
% Fo = F1, Fm = Fm'1

Fo = flstack(:,:,1);
Fm = flstack(:,:,2);
Fv = flstack(:,:,2) - flstack(:,:,1);

FvFm = Fv./Fm;
FvFm(FvFm < 0) = 0;



end

