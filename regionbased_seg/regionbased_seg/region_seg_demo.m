% Demo of "Region Based Active Contours"
%
% Example:
% seg_demo
%
% Coded by: Shawn Lankton (www.shawnlankton.com)

Img = imread('in000001.jpg');  %-- load the image
I = zeros(size(I,1),size(I,2));
m = zeros(size(I,1),size(I,2));          %-- create initial mask
m(200:400,200:630) = 1;
I = uint8(I);
I(200:400,200:630) = Img(200:400,200:630);
I = imresize(I,.5);  %-- make image smaller 
m = imresize(m,.5);  %     for fast computation

subplot(2,2,1); imshow(I); title('Input Image');
subplot(2,2,2); imshow(m); title('Initialization');
subplot(2,2,3); title('Segmentation');

seg = region_seg(I, m, 250); %-- Run segmentation

subplot(2,2,4); imshow(seg); title('Global Region-Based Segmentation');


