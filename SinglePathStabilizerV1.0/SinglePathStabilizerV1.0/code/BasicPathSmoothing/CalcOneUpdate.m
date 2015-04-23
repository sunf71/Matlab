function resultSet = CalcOneUpdate(updateSet,FList,weights,numFrames,num)

  resultSet = zeros(3,3,numFrames);
  HList = zeros(3,3,numFrames);
  HList(:,:,1) = eye(3,3);
  
  for i=1:numFrames-1
     HList(:,:,i+1) = inv(updateSet(:,:,i+1)) * FList(:,:,i+1) * updateSet(:,:,i);
  end
  
  postUpdateSet = CalcOneSmoothing(HList,weights,numFrames,num);
  
  for i=0:numFrames-1
     resultSet(:,:,i+1) = updateSet(:,:,i+1)*postUpdateSet(:,:,i+1);
  end
end