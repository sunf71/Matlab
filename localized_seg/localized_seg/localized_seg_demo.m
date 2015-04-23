% Demo of "Localized Region Based Active Contours"
% 
% Example:
% localized_seg_demo
%
% Coded by: Shawn Lankton (www.shawnlankton.com)

I = imread('monkey.png');         %-- load the image
m = false(size(I,1),size(I,2));   %-- create initial mask
m(37:213,89:227) = true;

I = imresize(I,.5);  %-- make image smaller 
m = imresize(m,.5);  %   for fast computation

subplot(2,2,1); imshow(I); title('Input Image');
subplot(2,2,2); imshow(m); title('Initialization');
subplot(2,2,3); title('Segmentation');

seg = localized_seg(I, m, 400);  %-- run segmentation

subplot(2,2,4); imshow(seg); title('Final Segmentation');




