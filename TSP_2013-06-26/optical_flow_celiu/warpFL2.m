function [warpI2,Imask]=warp(i2,vx,vy)
% warp i1 according to flow field in vx vy
[M,N,Z]=size(i2);
[x,y]=meshgrid(1:N,1:M);

warpI2 = zeros(size(i2));
for z=1:Z
    warpI2(:,:,z)=interp2(x,y,i2(:,:,z),x-vx,y-vy,'bicubic');
    %warpI2=interp2(x,y,i2,x+vx,y+vy,'linear');
    %warpI2=interp2(x,y,i2,x+vx,y+vy);
end

Imask = isnan(warpI2);
warpI2(Imask) = 0;