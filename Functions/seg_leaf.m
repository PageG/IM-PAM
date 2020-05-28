function [Lstack] = seg_leaf(image, F, n_leaves)
%seg_leaf Segment leaves in fluorescence image
%   Function to segment, label and extract leaves from within a
%   fluorescence image. 
%   If input image is a single n*m matrix, output is an 480*640*n array, where n = number of
%   leaves in the image.
%   If input is an n*m*p matrix, output is a 1*p cell array with each cell
%   containing an 480*640*n array, where n = number of
%   leaves in the image.

lsmin = 2000;
lbmin = 10;

A = F(:,:,1);
A(~isnan(A)) = 1;
A(isnan(A)) = 0;

[L,n] = bwlabel(A);

for i = 0:n%n-1
    tmp = L;
    tmp(tmp ~= i) = NaN;
    lsize(i+1) = nansum(nansum(tmp));
end

nleaves = nansum(lsize > lsmin);
idx = find(lsize > lsmin); idx1 = idx - 1;
dim = size(A);

% create simplifed label matrix
Lnew = L;
Lnew(~isnan(Lnew)) = NaN;
for i = 1:nleaves
    Lnew(L == idx1(i)) = L(L == idx1(i));
end

%% Find distances between pairs of leaf segments

binaryImage = im2bw(Lnew);

binaryImage = imfill(binaryImage, 'holes');

boundaries = bwboundaries(binaryImage);
numberOfBoundaries = size(boundaries, 1);

% exclude small boundaries
for i = 1:length(boundaries)
    sbound(i) = size(boundaries{i}, 1);
end
nbound = nansum(sbound > lbmin);
idx2 = find(sbound > lbmin); 
clear sbound

for k = 1 : nbound
    newBound{k} = boundaries{idx2(k)};
    sbound(k) = size(newBound{k}, 1);
end


% find the 'n_leaves' largest boundaries/objects
% calculate the distance of all objects to 'n_leaves' largest
% re-label to match closest large object

if nbound > n_leaves
    
    [~, I] = maxk(sbound, n_leaves);
    [~,II] = mink(sbound, length(sbound)-n_leaves);

    for i = 1:n_leaves
        pairs{i} = [repmat(I(i), length(II), 1), II']; % set up contrasts
    end

    % calculate triangular matrix of distances for n pairs
    % pairs = nchoosek(1:nbound,2);

    for j = 1:n_leaves
        for i = 1:length(pairs{j}(:,1))
            boundary1x = newBound{pairs{j}(i,1)}(:,2);
            boundary1y = newBound{pairs{j}(i,1)}(:,1);
            boundary2x = newBound{pairs{j}(i,2)}(:, 2);
            boundary2y = newBound{pairs{j}(i,2)}(:, 1);

            distx = (boundary1x' - boundary2x).^2;
            disty = (boundary1y' - boundary2y).^2;
            dista = sqrt(distx + disty);
            minDistance(i,j) = min(min(dista));

    %         indexofmin = find(dista == minDistance(i));
    %         [row, col] = ind2sub(size(dista),indexofmin(1)); % row/col of minimum value
    % 
    %         x(i,1) = boundary1x(col);
    %         y(i,1) = boundary1y(col);
    %         x(i,2) = boundary2x(row);
    %         y(i,2) = boundary2y(row);

    %         clear dista distx disty boundary1x boundary1y boundary2x boundary2y indexofmin row col
            clear dista distx disty boundary1x boundary1y boundary2x boundary2y
        end
    end

    % find the closest 'leaf' to each chunk
    for i = 1:length(II)
        [~,III(i)] = mink(minDistance(i,:), 1);
        nlabel(i,:) = [II(i) I(III(i))];
    end

    % I % idx of the predisignated leaves
    % II % idx of the unnasigned chunks
    % I(III) % idx of the matching leaf for chunk II


    % Re-label segments based on closest neighbor
    for i = 1:length(II)
        Lnew(Lnew == idx1(nlabel(i,1))) = idx1(nlabel(i,2));
    end
end

%% Re-label segments using top-to-bottom ordering

% matrix where columns = leaves, rows = rows
% values > 0 = leaf present in that row
if nbound > n_leaves
    for i = 1:n_leaves
        tmp = Lnew;
        tmp(Lnew ~= idx1(I(i))) = NaN;
        tmp1(:,i) = nansum(tmp');
        lorder(i, 1) = find(tmp1(:,i), 1, 'first');
        lorder(i, 2) = idx1(I(i));
    end
else
        for i = 1:nbound
        tmp = Lnew;
        tmp(Lnew ~= idx1(i)) = NaN;
        tmp1(:,i) = nansum(tmp');
        lorder(i, 1) = find(tmp1(:,i), 1, 'first');
        lorder(i, 2) = idx1(i);
        end
end

lorder = sortrows(lorder);
if nbound >= n_leaves
    lorder(:,3) = 1:n_leaves;
else
    lorder(:,3) = 1:nbound;
end

tmp = Lnew;

if nbound >= n_leaves
    for i = 1:n_leaves
        Lnew(tmp == lorder(i,2)) = lorder(i,3);
    end
else
    for i = 1:nbound
        Lnew(tmp == lorder(i,2)) = lorder(i,3);
    end
end
    
%% Figures ----------------------------------------------------------------

% % plot distance results
% subplot(1, 2, 1)
% imshow(binaryImage);
% axis on;
% hold on;
% for k = 1 : nbound
% 	thisBoundary = boundaries{idx2(k)};
%     newBound{k} = boundaries{idx2(k)};
% 	plot(thisBoundary(:,2), thisBoundary(:,1), 'r', 'LineWidth', 1);
% end
% hold off;
% 
% for i = 1:length(pairs(:,1))
%     if minDistance(i) < dist_crit
%         line(x(i,:), y(i,:), 'Color', 'g', 'Linewidth', 3);
%     else
%         line(x(i,:), y(i,:), 'Color', 'r', 'Linewidth', 3);
%     end
% end

% % plot relabelled results
% clear lsize
% for i = 0:n-length(mpairs)%n-1
%     tmp = Lnew;
%     tmp(tmp ~= i) = NaN;
%     lsize(i+1) = nansum(nansum(tmp));
% end
% 
% nleaves = nansum(lsize > lsmin);
% idx = find(lsize > lsmin); idx = idx - 1;
%% PLOTTING

% [nr,nc] = size(Lnew(:,:,1));
% % subplot(1, 2, 2)
% pcolor([Lnew nan(nr,1); nan(1,nc+1)]);
% shading flat;
% set(gca, 'ydir', 'reverse');
% % colormap hsv % jet, hsv

%%
% s = regionprops(binaryImage, 'Centroid');

% for i = 1:length(idx)
%     c = s(i).Centroid;
%     text(c(1), c(2), sprintf('%d', i), 'FontSize',14, ...
%         'HorizontalAlignment', 'center', ...
%         'VerticalAlignment', 'middle');
% end

%% -------------------------------------------------------------------------

tmp = size(image);


if nansum(~isnan(tmp)) < 3    
    Lstack = zeros(dim(1,1), dim(1,2), n_leaves);
    for i = 1:n_leaves
        A = image;
        A(Lnew ~= i) = NaN;
        Lstack(:,:,i) = A;
        clear A
    end
else
    Lstack = cell(1,tmp(1,3));
    for i = 1:tmp(1,3)        
        ltmp = zeros(dim(1,1), dim(1,2), n_leaves);
        for j = 1:n_leaves
            A = image(:,:,i);
            A(Lnew ~= j) = NaN;
            ltmp(:,:,j) = A;
            clear A
        end
        Lstack{1,i} = ltmp;
        clear ltmp
    end

end

