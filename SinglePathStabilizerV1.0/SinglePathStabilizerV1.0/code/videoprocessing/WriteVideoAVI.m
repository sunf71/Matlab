function WriteVideoAVI(filename,path,rate,frames)
myObj = VideoWriter(strcat(path,filename));
writerObj.FrameRate = rate;
open(myObj);
[~,~,~,num] = size(frames);
for i=1:num
    writeVideo(myObj,uint8(frames(:,:,:,i)));
end 
close(myObj);
 fprintf('\t Video OUTPUT DONE\t\n');
end