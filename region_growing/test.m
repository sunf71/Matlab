clear all;
I = im2double(imread('in000272.jpg'));
x = 99; y = 238;
J = regiongrowing(I,x,y,0.1);
figure,imshow(I+J);