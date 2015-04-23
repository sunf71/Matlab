function showIM=video_horizontal(im1,im2,gap)
    [row1,col1] = size(im1);
    [row2,col2] = size(im2);
    
    if row1>row2
         im2(row2+1:row1,:) = 0;
         gappart=zeros(row1,gap,3);
         im = [im1 gappart im2]; 
    end
    if row2>row1
         im1(row1+1:row2,:) = 0;
         gappart=zeros(row2,gap,3);
         im = [im1 gappart im2];
    end
    if row2==row1
         gappart=zeros(row1,gap,3);
         im = [im1 gappart im2];   
    end
    
    showIM=im;
    
end

