function data = SP_img2data(img, w)

xdim = size(img,1);
ydim = size(img,2);
num_pix = (2*w+xdim)*(2*w+ydim);
data = zeros(5, num_pix);
[yg xg] = meshgrid(-w:ydim-1+w, -w:xdim-1+w);
data(1,:) = xg(:)';
data(2,:) = yg(:)';

% lab = -ones(size(xg,1), size(xg,2), 3);
% lab(w+1:w+xdim, w+1:w+ydim,:) = rgb2lab(img);

if (~isa(img, 'uint8'))
    img = im2uint8(img);
end
lab = double(rgb2lab(img));
lab = holdpad(lab, xdim+2*w, ydim+2*w);

data(3,:) = reshape(lab(:,:,1),1,num_pix);
data(4,:) = reshape(lab(:,:,2),1,num_pix);
data(5,:) = reshape(lab(:,:,3),1,num_pix);


% cform = makecform('srgb2lab');
% img_lab = lab2double(applycform(img,cform));
% data(3,:) = reshape(img_lab(:,:,1),1,num_pix)*2.55;
% data(4,:) = reshape(img_lab(:,:,2),1,num_pix)+ 128;
% data(5,:) = reshape(img_lab(:,:,3),1,num_pix)+ 128;



