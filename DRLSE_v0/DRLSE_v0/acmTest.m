%  This Matlab code demonstrates an edge-based active contour model as an application of 
%  the Distance Regularized Level Set Evolution (DRLSE) formulation in the following paper:
%
%  C. Li, C. Xu, C. Gui, M. D. Fox, "Distance Regularized Level Set Evolution and Its Application to Image Segmentation", 
%     IEEE Trans. Image Processing, vol. 19 (12), pp. 3243-3254, 2010.
%
% Author: Chunming Li, all rights reserved
% E-mail: lchunming@gmail.com   
%         li_chunming@hotmail.com 
% URL:  http://www.imagecomputing.org/~cmli//

clear all;
close all;
id = 19;
imgName = sprintf('in%06d.jpg',id);
Img=imread(imgName);
Img = rgb2gray(Img);
Img = imresize(Img,.5);
Img=double(Img(:,:,1));
imgName = sprintf('bin%06d.png',id);
InitContour = imread(imgName);
InitContour = imresize(InitContour,.5);
InitContour = rgb2gray(InitContour);
mask = false(size(InitContour));
mask(InitContour==255) = true;

sigma=1.5;   % scale parameter in Gaussian kernel for smoothing.
epsilon=1.5; % the papramater in the definition of smoothed Dirac function
mu=0.04;     % coefficient of the internal (penalizing) energy term P(\phi)
lambda=5;    % coefficient of the weighted length term Lg(\phi)
alf=3;     % coefficient of the weighted area term Ag(\phi);
N=100;%number of iterations
c0=4;        % the constant value used to define binary level set function;
PlotRate=10;

u = LevelSetEvolutionWithoutReinitialization(Img,sigma,epsilon,mu,lambda,alf,c0,N,PlotRate,mask);
contourP = contourf(u,[0 0]);
contourP = uint16(contourP+0.5);
result = zeros(size(Img));
result(contourP) = 255;
result = uint8(result);
fileName = sprintf('bin%06d.jpg',id);
imwrite(result,fileName);
% %% parameter setting
% timestep=1;  % time step
% mu=0.2/timestep;  % coefficient of the distance regularization term R(phi)
% iter_inner=5;
% iter_outer=20;
% lambda=5; % coefficient of the weighted length term L(phi)
% alfa=-3;  % coefficient of the weighted area term A(phi)
% epsilon=1.5; % papramater that specifies the width of the DiracDelta function
% 
% sigma=.8;    % scale parameter in Gaussian kernel
% G=fspecial('gaussian',15,sigma); % Caussian kernel
% Img_smooth=conv2(Img,G,'same');  % smooth image by Gaussiin convolution
% [Ix,Iy]=gradient(Img_smooth);
% f=Ix.^2+Iy.^2;
% g=1./(1+f);  % edge indicator function.
% 
% % initialize LSF as binary step function
% c0=2;
% initialLSF = c0*ones(size(Img));
% % generate the initial region R0 as two rectangles
% 
% initialLSF(InitContour==255) = -c0;
% phi=initialLSF;
% 
% figure(1);
% mesh(-phi);   % for a better view, the LSF is displayed upside down
% hold on;  contour(phi, [0,0], 'r','LineWidth',2);
% title('Initial level set function');
% view([-80 35]);
% 
% figure(2);
% imagesc(Img,[0, 255]); axis off; axis equal; colormap(gray); hold on;  contour(phi, [0,0], 'r');
% title('Initial zero level contour');
% pause(0.5);
% 
% potential=2;  
% if potential ==1
%     potentialFunction = 'single-well';  % use single well potential p1(s)=0.5*(s-1)^2, which is good for region-based model 
% elseif potential == 2
%     potentialFunction = 'double-well';  % use double-well potential in Eq. (16), which is good for both edge and region based models
% else
%     potentialFunction = 'double-well';  % default choice of potential function
% end  
% 
% % start level set evolution
% for n=1:iter_outer
%     phi = drlse_edge(phi, g, lambda, mu, alfa, epsilon, timestep, iter_inner, potentialFunction);    
%     if mod(n,2)==0
%         figure(2);
%         imagesc(Img,[0, 255]); axis off; axis equal; colormap(gray); hold on;  contour(phi, [0,0], 'r');
%     end
% end
% 
% % refine the zero level contour by further level set evolution with alfa=0
% alfa=0;
% iter_refine = 10;
% phi = drlse_edge(phi, g, lambda, mu, alfa, epsilon, timestep, iter_inner, potentialFunction);
% 
% finalLSF=phi;
% figure(2);
% imagesc(Img,[0, 255]); axis off; axis equal; colormap(gray); hold on;  contour(phi, [0,0], 'r');
% hold on;  contour(phi, [0,0], 'r');
% str=['Final zero level contour, ', num2str(iter_outer*iter_inner+iter_refine), ' iterations'];
% title(str);
% 
% figure;
% mesh(-finalLSF); % for a better view, the LSF is displayed upside down
% hold on;  contour(phi, [0,0], 'r','LineWidth',2);
% view([-80 35]);
% str=['Final level set function, ', num2str(iter_outer*iter_inner+iter_refine), ' iterations'];
% title(str);
% axis on;
% [nrow, ncol]=size(Img);
% axis([1 ncol 1 nrow -5 5]);
% set(gca,'ZTick',[-3:1:3]);
% set(gca,'FontSize',14)

