function updateSet = CalcOneSmoothing(motionSet,weights,numFrames,num)

    updateSet = zeros(3,3,numFrames);
    
    motionSet_Inv = zeros(3,3,numFrames);
    for i=0:numFrames-1
        motionSet_Inv(:,:,i+1) = inv(motionSet(:,:,i+1));
    end
    
    length = num*2+1;
    windowMotions = zeros(3,3,length);
    nbHood = zeros(length,1);

    
    
for t=0:numFrames-1
        for k=0:length-1
            windowMotions(:,:,k+1) = eye(3);
            nbHood(k+1) = 0;
        end
        
        k = num-1;
        for curIdx=t:-1:t-num+1
            if curIdx>0
                windowMotions(:,:,k+1) = windowMotions(:,:,k+1+1) * motionSet(:,:,curIdx+1);
                nbHood(k+1) = 1;
            end
            k = k-1;
        end
        
        k = num+1;
        for curIdx=t+1:t+num
            if curIdx<numFrames
                windowMotions(:,:,k+1) = windowMotions(:,:,k-1+1)*motionSet_Inv(:,:,curIdx+1);
                nbHood(k+1) = 1;
            end
            k = k+1;
        end
        
        sumMotion = zeros(3);
        sumWeight = 0;
        for i=0:length-1
           if nbHood(i+1)
              sumMotion = sumMotion +  weights(i+1) .* windowMotions(:,:,i+1);
              sumWeight = sumWeight + weights(i+1);
           end
        end
        
        updateSet(:,:,t+1) = sumMotion ./ sumWeight;
end
    
end

