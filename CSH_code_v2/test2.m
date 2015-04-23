close all
clear all
clc

AddPaths
% preparing the images
img1 = 'in000009.jpg';
img2 = 'in000010.jpg';

A = imread(img1);
B = imread(img2);
grayA = rgb2gray(A);
grayB = rgb2gray(B);
[hA wA dA] = size(A);
[hB wB dB] = size(B);

mpa = floor(hB*wB/1000) / 1000;
mpb = floor(hB*wB/1000) / 1000;

MP_A_Str = [num2str(mpa) ' MP'];
MP_B_Str = [num2str(mpb) ' MP'];

fprintf('Image A: %s, size = %s\r\n', img1 , MP_A_Str);
fprintf('Image B: %s, size = %s\r\n', img2 , MP_B_Str);
fprintf('1. Runing CSH_nn example with default parameter values\r\n');
width = 16;
tic;
%%%%%% CSH RUN %%%%%%%

CSH_ann = CSH_nn(A,B,width,5);
CSH_TOC = toc
fid=fopen('data.txt','w+');
for i=1 : 500  
    x = floor(rand*wA)+1;
    y = floor(rand*hA)+1;
    wx = CSH_ann(y,x,1)+1;
    wy =  CSH_ann(y,x,2)+1;
    diff = abs(grayA(y,x,1) - grayB(wy,wx,1));
    if (diff <20)
        fprintf(fid,'%d %d %d %d \r\n',x,y,wx,wy);
    end
    
end
fclose(fid);