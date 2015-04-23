function FList = CalcMotion(name,input_frames)

[h,w,c,numFrames] = size(input_frames(:,:,:,:));

FList = zeros(3,3,numFrames);
FList(:,:,1) = eye(3);

    for i=2:numFrames
         s = sprintf('compute homography %d / %d \n',i-1,numFrames);
         s
        FList(:,:,i) = EstimateHomography_SURF_RANSAC(input_frames(:,:,:,i),input_frames(:,:,:,i-1));
    end
    save(strcat(name,'\homography'),'FList');
end