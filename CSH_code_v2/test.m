close all
clear all
clc

AddPaths

fprintf('CSH algorithm example script!!!\r\n');
fprintf('*******************************\r\n');


% preparing the images
img1 = 'in000289.jpg';
img2 = 'in000290.jpg';

A = imread(img1);
B = imread(img2);

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
CSH_TIC = tic;
%%%%%% CSH RUN %%%%%%%
CSH_ann = CSH_nn(A,B,width,5);
%%%%%%%%%%%%%%%%%%%%%%
CSH_TOC = toc(CSH_TIC);
fprintf('   CSH_nn elapsed time: %.3f[s]\r\n' , CSH_TOC);
homography = [0.995,-0.00177,0.73323;0.0009,0.9944785,0.322796;4.1425e-6,-1.502e-5,1];
for i=1:5
%CSH_ann(mod(rand(),hA),mod(rand(),wA))
mod(rand(),hA)
end
mask = zeros(hA,wA);
dAnn = double(CSH_ann);
for i=1:hA
    for j=1:wA
        x = homography(1,1)*j + homography(1,2)*i + homography(1,3);
		y = homography(2,1)*j + homography(2,2)*i + homography(2,3);
		w = homography(3,1)*j + homography(3,2)*i + homography(3,3);
		x = x/w;
		y =  y/w;
        if abs(x-dAnn(i,j,1)) + abs(y-dAnn(i,j,2)) > 1
        mask(i,j) = 255;
        end
    end
end
imshow(mask);

patchA = zeros(width*width*dA,hA*wA);
patchB = zeros(width*width*dB,hA*wA);
for i=1:hA
    for j=1:wA
        pos(i,j,1) = j;
        pos(i,j,2) = i;
        if i+width < hA && j+width < wA
        for y=1:width-1
            for x=1:width-1
                for k=1:3
                    patchA((x+y*width)*3+k,(i*wA + j)) = A(i+y,j+x,k);
                    patchB((x+y*width)*3+k,(i*wA + j)) = B(i+y,j+x,k);
                end
            end
        end
        end
    end 
end

patchParms = struct('hA',hA,'wA',wA,'dA',dA,'hB',hB,'wB',wB,'dB',dB);
 %PlotExampleResults(A,B,CSH_Mapping,width,K_of_KNN,bMask,experimentName,patch_mode,A_patch,B_patch,patch_params,interactive,...
  %                                                          compensateOrientations,winnerOrientations_A,winnerOrientations_B,A_patch_normal,B_patch_normal)
PlotExampleResults(A,B,CSH_ann,width,1,[],'default CSH',1,patchA,patchB,patchParms,1);