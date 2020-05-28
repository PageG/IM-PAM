function [Fo, Fm, Fv, FvFm, Fmp, Fop, F, YII, NPQ, qN, qP, qL, YNPQ, YNO, ETR] = im_pam_tiff(file,PAR)
%im_pam_tiff Import TIFF stack from Waltz Imaging PAM
%   Use this function to import and process and TIFF stack exported from
%   ImagingWin software (tested with ImagingWin v2.46i, with data collected using an M-series Imaging PAM
%   fitted with an IMAG-MIN/B measuring head, Walz GmbH, Effeltrich, Germany. Outputs are calculated from
%   formalae described in the Imaging PAM manual.
%   PAR can either be a single value (as for a dark-light induction curve)
%   or a vector (as for a rapid light curve).

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

for i = 1:((num_images-4)/2)
    F(:,:,i) = flstack(:,:,((i*2)+3));
end

for i = 1:((num_images-4)/2)
    Fmp(:,:,i) = flstack(:,:,((i*2)+4));
end

% Fo' = Fo/ (Fv/Fm + Fo/Fm')
for i = 1:((num_images-4)/2)
    Fop(:,:,i) = Fo ./ ((Fv ./ Fm) + (Fo ./ Fmp(:,:,i)));
end

% Y(II) = (Fm'-F)/Fm'
for i = 1:((num_images-4)/2)
    A = (Fmp(:,:,i) - F(:,:,i)) ./ Fmp(:,:,i);
    A(A < 0) = NaN;
    A(A > 1) = NaN;
    YII(:,:,i) = A;
    clear A
end

% NPQ = (Fm-Fm')/Fm'
for i = 1:((num_images-4)/2)
    NPQ(:,:,i) = (Fm - Fmp(:,:,i)) ./ Fmp(:,:,i);
end

% qN = (Fm-Fm')/(Fm-Fo')
for i = 1:((num_images-4)/2)
    qN(:,:,i) = abs(Fm - Fmp(:,:,i)) ./ (Fm - Fop(:,:,i));
    A = qN(:,:,i);
    A(A < 0) = NaN;
    A(A > 1) = NaN;
    qN(:,:,i) = A;
    clear A
end

% qP = (Fm'-F)/(Fm'-Fo')
for i = 1:((num_images-4)/2)
    qP(:,:,i) = (Fmp(:,:,i) - F(:,:,i)) ./ (Fmp(:,:,i) - Fop(:,:,i));
    A = qP(:,:,i);
    A(A < 0) = NaN;
    A(A > 1) = NaN;
    qP(:,:,i) = A;
    clear A
end

% qL = (Fm'-F)/(Fm'-Fo') x Fo'/F = qP x Fo'/F
for i = 1:((num_images-4)/2)
    qL(:,:,i) = qP(:,:,i) .* (Fop(:,:,i) ./ F(:,:,i));
end

% Y(NPQ) = 1 - Y(II) - 1/(NPQ+1+qL(Fm/Fo-1))
for i = 1:((num_images-4)/2)
    A = 1 - YII(:,:,i) - 1 ./(NPQ(:,:,i) + 1 + qL(:,:,i) .* (Fm ./ Fo - 1));
    A(A < 0) = NaN;
    YNPQ(:,:,i) = A;
    clear A
end

% Y(NO) = 1/(NPQ+1+qL(Fm/Fo-1))
for i = 1:((num_images-4)/2)
    A = 1 ./(NPQ(:,:,i) + 1 + qL(:,:,i) .* (Fm ./ Fo - 1));
    A(A < 0) = NaN;
    YNO(:,:,i) = A;
    clear A
end

% ETR = 0.5 x Yield x PAR x 0.84
for i = 1:((num_images-4)/2)
    if nansum(~isnan(PAR)) > 1
        ETR(:,:,i) = 0.5 .* YII(:,:,i) .* PAR(1,i) .* 0.84;
    else
        ETR(:,:,i) = 0.5 .* YII(:,:,i) .* PAR .* 0.84;
    end

end


