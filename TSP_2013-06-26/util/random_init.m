function [z, num_z] = random_init(xdim, ydim, w, K)

rxdim = xdim - 2*w;
rydim = ydim - 2*w;

N = rxdim*rydim;

%centers = randperm(N, K);
centers = randsample(N, K);

[centers_x, centers_y] = ind2sub([rxdim,rydim], centers);
centers_x = centers_x + w;
centers_y = centers_y + w;
centers = sub2ind([xdim,ydim], centers_x, centers_y);

z = int32(random_initIMPORT(xdim,ydim, centers));

% [centers_x, centers_y] = ind2sub([xdim, ydim], centers);
% 
% [yg xg] = meshgrid(1:ydim, 1:xdim);
% 
% distances = zeros(N,K);
% for k=1:K
%     distances(:,k) = (yg(:)-centers_y(k)).^2 + (xg(:)-centers_x(k)).^2;
% end
% 
% [~, z] = min(distances, [], 2);
% z = uint32(reshape(z, xdim, ydim));
% z = Util_labelconnect(z);
% 
num_z = double(max(z(:)));