clear;
img=imread('in000010.jpg');
lab=rgb2lab(img); % You need some function RGB2Lab to perform the conversion
l = lab(:,:,1);
a = lab(:,:,2);
b = lab(:,:,3);
rr= img(:,:,1);
gg = img(:,:,2);
bb = img(:,:,3);
%dataCell = {l,a,b};
dataCell = {rr,gg,bb};
[salMat,salMatInd]=saliencyMeasure(dataCell);