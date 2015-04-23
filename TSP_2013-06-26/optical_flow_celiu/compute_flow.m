function [angle, vx, vy] = compute_flow(im1, im2)

if (max(im1(:))>1)
    im1 = double(im1) / 255;
end
if (max(im2(:))>1)
    im2 = double(im2) / 255;
end

% set optical flow parameters (see Coarse2FineTwoFrames.m for the definition of the parameters)
alpha = 0.012;
%alpha = 0.05;
ratio = 0.75;
%ratio = 1;
minWidth = 50;
nOuterFPIterations = 7;
nInnerFPIterations = 1;
nSORIterations = 30;

para = [alpha,ratio,minWidth,nOuterFPIterations,nInnerFPIterations,nSORIterations];

% this is the core part of calling the mexed dll file for computing optical flow
% it also returns the time that is needed for two-frame estimation
[vx,vy,warpI2] = Coarse2FineTwoFrames(im1,im2,para);
angle = atan2(vy, vx);