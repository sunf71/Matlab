clear
close all
f = 'a';
ext = 'JPG';
img1 = imread(['hall1.' ext]);
img2 = imread(['hall2.' ext]);
% img3 = imread('b3.jpg');

imMosaic(img2,img1,1);
% img0 = imMosaic(img1,img0,1);