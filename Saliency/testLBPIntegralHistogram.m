clear;
img = imread('in000010.jpg');
gray = rgb2gray(img);
mapping = getmapping(8,'u2');
lbpImg = lbp(gray,1,8,mapping,'');
blbpImg = zeros(size(lbpImg,1)+2,size(lbpImg,2)+2);
blbpImg(1,2:size(lbpImg,2)+1)=lbpImg(1,:);
blbpImg(2:size(lbpImg,1)+1,1)=lbpImg(:,1);
blbpImg(size(lbpImg,1)+2,2:size(lbpImg,2)+1)=lbpImg(size(lbpImg,1),:);
blbpImg(2:size(lbpImg,1)+1,size(lbpImg,2)+2)=lbpImg(:,size(lbpImg,2));
blbpImg(1,size(blbpImg,2)) = lbpImg(1,1);
blbpImg(1,1)  = lbpImg(1,1);
blbpImg(1,size(blbpImg,2)) = lbpImg(1,1);
blbpImg(size(blbpImg,1),1)  = lbpImg(size(lbpImg,1),1);
blbpImg(size(blbpImg,1),size(blbpImg,2))  = lbpImg(size(lbpImg,1),size(lbpImg,2));
blbpImg(2:size(blbpImg,1)-1,2:size(blbpImg,2)-1) = lbpImg;

minLBP = min(min(blbpImg));
maxLBP = max(max(blbpImg));
height = size(img,1);
width = size(img,2);

[intHist,intInd]=formLBPIntegralHistogram(lbpImg);

h11 = 10;h12 = 10; h21 = 11;h22 = 11;

 hist1 = intHist(h21,h22,:) + intHist(h11-1,h12-1,:) - intHist(h11-1,h22,:) - intHist(h21,h12-1,:);
 hist1 = hist1(:)';
 window = gray(h11-1:h21+1,h12-1:h22+1);
 hist2 = lbp(window,1,8,mapping,'h');
 fprintf('total %d pixels in hist2, and %d pixels in hist1\n',sum(hist2),sum(hist1));
 subplot(2,1,1), stem(hist1);
 subplot(2,1,2), stem(hist2);