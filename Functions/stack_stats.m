function [Smean, Ssd] = stack_stats(stack)
%   stack_stats Compute stats on a fluorescence image stack
%   Input is a segmented image stack (see 'seg_leaf.m')
%   Output is mean and standard error for each leaf in each image

tmp = stack{1,1};

for j = 1:length(stack)
    for i = 1:size(tmp,3)
       Smean(i,j) = nanmean(nanmean(stack{1,j}(:,:,i)));
       Ssd(i,j) = std(reshape(stack{1,j}(:,:,i), 1, 640*480), 'omitna');

    end
end

% figure()
% for i = 1:size(Smean,1)
%     plot(Smean(i,:))
%     hold on
% end

end

