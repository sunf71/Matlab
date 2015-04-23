function [padded] = holdpad(x, dimx, dimy)

C = size(x,3);

padded = zeros(dimx,dimy,C, class(x));

xstart = ceil((dimx-size(x,1))/2)+1;
xstop = xstart+size(x,1)-1;
ystart = ceil((dimy-size(x,2))/2)+1;
ystop = ystart+size(x,2)-1;

for c=1:C
    padded(xstart:xstop,ystart:ystop,c) = x(:,:,c);
    padded(xstart:xstop, 1:ystart-1,c) = padded(xstart:xstop, ystart,c) * ones(1, ystart-1);
    padded(xstart:xstop, ystop+1:size(padded,2),c) = padded(xstart:xstop, ystop,c) * ones(1, size(padded,2)-ystop);
    padded(1:xstart-1, :,c) = ones(xstart-1, 1) * padded(xstart, :,c);
    padded(xstop+1:size(padded,1), :,c) = ones(size(padded,1)-xstop, 1) * padded(xstop, :,c);
end