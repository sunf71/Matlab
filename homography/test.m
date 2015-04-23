close all
clear all
clc

% preparing the images
img1 = 'in000009.jpg';
img2 = 'in000010.jpg';

A = imread(img1);
B = imread(img2);
grayA = rgb2gray(A);
grayB = rgb2gray(B);

[ h C ] = homography( grayA, grayB);
%imshow('warped',C);
