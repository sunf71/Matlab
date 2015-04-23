function warpI = HomographyWarp(I,H )
    [h,w,~] = size(I);
    szIm = [h,w];
    [x,y] = meshgrid(1:w,1:h);
    pix = [x(:)'; y(:)'];
    hPixels = [ pix; ones(1,prod(szIm))];
    hScene  = H*hPixels;
    xprime=(hScene(1,:)./(hScene(3,:)))';
    yprime=(hScene(2,:)./(hScene(3,:)))';

    xprime = reshape(xprime, szIm);
    yprime = reshape(yprime, szIm);
    
    warpI(:,:,1) = interp2(I(:,:,1),xprime,yprime,'cubic');
    warpI(:,:,2) = interp2(I(:,:,2),xprime,yprime,'cubic');
    warpI(:,:,3) = interp2(I(:,:,3),xprime,yprime,'cubic');
end

