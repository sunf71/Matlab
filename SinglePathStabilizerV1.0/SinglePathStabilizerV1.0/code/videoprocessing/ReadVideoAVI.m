function [frameNum,frameRate,frames]=ReadVideoAVI(fileName,path)

name = strcat(path,fileName);
obj = VideoReader(name);
frameNum = obj.NumberOfFrames;
frameRate = obj.FrameRate;
frameNum=frameNum-1;  
frames = read(obj,[1,frameNum]);


end