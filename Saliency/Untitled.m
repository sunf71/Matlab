clear;
img = imread('in000005.jpg');
labImg = rgb2lab(img);
gray = rgb2gray(img);
mapping = getmapping(8,'u2');
lbpImg = lbp(gray,1,8,mapping,'');
% blbpImg = zeros(size(lbpImg,1)+2,size(lbpImg,2)+2);
% blbpImg(1,2:size(lbpImg,2)+1)=lbpImg(1,:);
% blbpImg(2:size(lbpImg,1)+1,1)=lbpImg(:,1);
% blbpImg(size(lbpImg,1)+2,2:size(lbpImg,2)+1)=lbpImg(size(lbpImg,1),:);
% blbpImg(2:size(lbpImg,1)+1,size(lbpImg,2)+2)=lbpImg(:,size(lbpImg,2));
% blbpImg(1,size(blbpImg,2)) = lbpImg(1,1);
% blbpImg(1,1)  = lbpImg(1,1);
% blbpImg(1,size(blbpImg,2)) = lbpImg(1,1);
% blbpImg(size(blbpImg,1),1)  = lbpImg(size(lbpImg,1),1);
% blbpImg(size(blbpImg,1),size(blbpImg,2))  = lbpImg(size(lbpImg,1),size(lbpImg,2));
% blbpImg(2:size(blbpImg,1)-1,2:size(blbpImg,2)-1) = lbpImg;
% minLBP = min(min(blbpImg));
% maxLBP = max(max(blbpImg));
height = size(img,1);
width = size(img,2);

[intHist,intInd]=formLBPIntegralHistogram(lbpImg);
% [intHist{2},dataCell{2}]=formIntegralHistogram(lab(:,:,2),60,-128,128);
% [intHist{3},dataCell{3}]=formIntegralHistogram(lab(:,:,3),60,-128,128);
step = floor(0.04*width);
sal = zeros(height,width);
ssdConfidence = 0.6;
% wh = 200;
% ww = 79;
 wh = 150;
 ww = 80;
windowRows=[0.25 0.3 0.5 0.7]; % Window row size relative to image size max(imRow,imCol))
windowCols=[0.1 0.3 0.5 0.4]; % Window column size relative to image size max(imRow,imCol))
sampleStep=[0.01 0.015 0.03 0.04]; % Window sampling step relative to image size max(imRow,imCol)) (one for each window size) 
maxH = max(width,height);
whs = floor(windowRows*maxH);
wws = floor(windowCols*maxH);
steps = floor(sampleStep*max(whs,wws));

minWindow = zeros(wh+1,ww+1);
tic;
% sh = 1;
% sw = 1;
sh = 150;
sw = 280;
maxSal = 0;
maxWindows =zeros(wh+1,ww+1);
maxWindowt =zeros(wh+1,ww+1);

for ii=1:step:480
    for jj = 1:step:640
        if ii + wh < height && jj + ww <width
            ll = [ii jj];
            srect=[ll(2),ll(1),ll(2)+ww,ll(1)+wh];
            window1 = double(img(ll(1):ll(1)+wh,ll(2):ll(2)+ww,:));
            lab1 = labImg(ll(1):ll(1)+wh,ll(2):ll(2)+ww,:);
          
            %imshow(uint8(window1));
            
            hist1 = intHist(ll(1)+wh,ll(2)+ww,:) + intHist(ii,jj,:) - intHist(ii,jj+ww,:) - intHist(ii+wh,jj,:);
            hist1 = double(hist1(:)')/sum(hist1);
            minSSD = 1e10;
            minDiff = zeros(wh,ww);
            c = 0;
            for i=1:step:height
                for j=1:step:width
                    if j+ww< width && i +wh <height && (j>srect(3) || i > srect(4) || srect(1) > j+ww || srect(2) > i+wh)
                        window  = double(img(i:i+wh,j:j+ww,:));
                        lab = labImg(i:i+wh,j:j+ww,:);
                        c=c+1;
                        hist2 = intHist(i+wh,j+ww,:) + intHist(i,j,:) - intHist(i,j+ww,:) - intHist(i+wh,j,:);
                        hist2 = double(hist2(:)')/sum(hist2);
                        diff = abs(lab1-lab);
                        diff = (diff(:,:,1)/100 + diff(:,:,2)/255 + diff(:,:,3)/255)/3;
                        ssd = sum(sum(diff))/wh/ww;
                        %ssd = pdist2(hist1(:)',hist2(:)','emd');
                        histDist = bhattacharyya(hist1,hist2);
                        patchDist = ssd*ssdConfidence + histDist*(1-ssdConfidence);
                        if patchDist < minSSD
                            minSSD = patchDist;
                            minWindow = window;
%                             minDiff = diff;
                        end
                        
                    end
                end
            end
            if (minSSD > maxSal)
                maxSal = minSSD;
                maxWindows = window1;
                maxWindowt = minWindow;
            end
            sal(ll(1):ll(1)+wh,ll(2):ll(2)+ww,:) = max(sal(ll(1):ll(1)+wh,ll(2):ll(2)+ww,:),minSSD);
        end
    end
end

toc
fprintf('%d window\n',c);

fprintf('minSSD = %f\n',minSSD);
min = min(min(sal));
max = max(max(sal));
nsal = (sal-min)/(max-min);

% nsal = im2bw(nsal,0.8);
nsal = nsal*255;
% subplot(1,3,1);
imshow(uint8(nsal));
% subplot(1,3,2);
% imshow(uint8(window1));
% subplot(1,3,3);
% imshow(uint8(minWindow));
% imwrite(uint8(maxWindows),'s.jpg');
% imwrite(uint8(maxWindowt),'t.jpg');
%imshow(window2);
