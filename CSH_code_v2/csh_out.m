%csh output
istart = 1;
iend = 19;
preImg = [];
for i=istart : iend
    imgName = sprintf('in%06d.jpg',i);
    img = imread(imgName);
    %imgName = sprintf('in%06d.jpg',i-1);
    if size(preImg,1) == 0
        preImg = img;
    end
    ann = CSH_nn(img,preImg);
    dataName = sprintf('flow%06d.txt',i);
    dlmwrite(dataName,ann);
    fid=fopen(dataName,'w+');
    height = size(img,1);
    width = size(img,2);   
    for m=1:height
        for n = 1:width-1            
            fprintf(fid,'%d,%d,',ann(m,n,1),ann(m,n,2));
        end
        fprintf(fid,'%d,%d\n',ann(m,width,1),ann(m,width,2));
    end
    fclose(fid);
%     height = size(img,1);
%     width = size(img,2);
%     dst = uint8(zeros(height,width,3));
%     for m=1:height
%         for n = 1:width
%             dst(m,n,:) = preImg(ann(m,n,2),ann(m,n,1),:);
%         end
%     end
%     imshow(dst);
   preImg = img;
end