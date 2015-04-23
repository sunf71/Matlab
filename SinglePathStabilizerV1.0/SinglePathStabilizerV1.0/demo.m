clear all;
clc;

addpath('code/MotionEstimation/');
addpath('code/MotionEstimation/RANSAC');
addpath('code/BasicPathSmoothing/');
addpath('code/videoprocessing');
folder = 'data/1/';            
InputVideoName='input.avi';


[numFrames,frameRate,frames]=ReadVideoAVI(InputVideoName,folder);

% % Alternatively, load video frames as image sequence
% Files = dir(strcat(folder,'frames/*.jpg'));
% numFrames = length(Files);
% for i = 1:numFrames;
%    frames(:,:,:,i) = im2double(imread(strcat(folder,'frames/',Files(i).name)));
% end
% frameRate = 30; % default frame rate


%estimate homography between neighboring frames
[h,w,c,~] = size(frames(:,:,:,:));

if exist(strcat(folder,'/homography.mat'),'file')==0
    CalcMotion(folder,frames);
end
FList = struct2array(load(strcat(folder,'/homography.mat')));

%single path smoothing
disp('path smoothing...');
windowsize = 15;
weights = CalcNeighborWeight(windowsize,sqrt(windowsize));
updateSet = zeros(3,3,numFrames);
for i=1:numFrames
   updateSet(:,:,i) = eye(3);
end

iterations = 20;
for i=1:iterations
   updateSet = CalcOneUpdate(updateSet,FList,weights,numFrames,windowsize);
end
disp('[DONE]');

%output stabilized frames
disp('rendering...');
for i=1:numFrames
  fprintf('%3d / %3d\n',i,numFrames);
  result(:,:,:,i) = HomographyWarp(frames(:,:,:,i),inv(updateSet(:,:,i)));
  ComparisionIM(:,:,:,i) = video_horizontal(frames(:,:,:,i),result(:,:,:,i),20);
end

WriteVideoAVI('output.avi',folder,frameRate,result);
WriteVideoAVI('comparison.avi',folder,frameRate,ComparisionIM);

disp('[DONE]');












