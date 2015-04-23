function motion = CalcMotionFor2Frame(motionSet,srcIdx,dstIdx)
    leftIdx = min(srcIdx,dstIdx);
    rightIdx = max(srcIdx,dstIdx);
    
    tmpMotion = eye(3);
    for t = rightIdx:-1:leftIdx+1
        tmpMotion = tmpMotion*motionSet(t+1);
    end
    
    if srcIdx >= dstIdx
        motion = tmpMotion;
    else
        motion = inv(tmpMotion);
    end
end