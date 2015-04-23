function [z, num_z] = random_prop(IMG)

xdim = IMG.xdim;
ydim = IMG.ydim;

N = xdim*ydim;
K = IMG.K;

centers_x = zeros(1,K);
centers_y = zeros(1,K);
for i=1:K
    centers = IMG.SP(i).p_mu + IMG.SP(i).v;
    centers_x(i) = centers(1);
    centers_y(i) = centers(2);
end

[yg xg] = meshgrid(1:ydim, 1:xdim);

distances = zeros(N,K);
for k=1:K
    distances(:,k) = (yg(:)-centers_y(k)).^2 + (xg(:)-centers_x(k)).^2;
end

[~, z] = min(distances, [], 2);
z = uint32(reshape(z, xdim, ydim));

for k=1:K
    [cc, Ncc] = bwlabel(z==k,4);
    if (Ncc>1)
        num_pixels = hist(cc(:),0:Ncc);
        num_pixels = num_pixels(2:end);
        [~, indices] = sort(num_pixels, 'descend');
        for c=2:Ncc
            [xs,ys] = find(cc==indices(c));
            for i=1:numel(xs)
                x = xs(i);
                y = ys(i);
                if (x>1 && z(x-1,y)~=k)
                    z(x,y) = z(x-1,y);
                elseif (y>1 && z(x,y-1)~=k)
                    z(x,y) = z(x,y-1);
                elseif (x<xdim && z(x+1,y)~=k)
                    z(x,y) = z(x+1,y);
                else %if (y<ydim)
                    z(x,y) = z(x,y+1);
                end
            end
        end
    end
end


num_z = double(max(K,max(z(:))));
