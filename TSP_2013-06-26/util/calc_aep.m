function [aep] = calc_aep(g, gt, mask)
ep = sqrt(sum((g - gt).^2,3));

if (~exist('mask','var') || isempty(mask))
    aep = mean(ep(:));
else
    aep = mean(ep(mask));
end