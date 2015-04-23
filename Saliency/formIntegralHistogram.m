%% Make integral histogram image
function [intHist,binInd]=formIntegralHistogram(dataMat,numBin,minVal,maxVal)
% currently only 1 histograms

% Compute bin indices
%binInd=ceil(numBin*(double(dataMat)+1e-5)/maxVal);
%binInd(binInd>numBin)=numBin;

if isempty(dataMat)
    error('Datamatrix is empty');
end

% If only one numBin/minVal/maxVal is given, but dataMat has more feature, use same numBin/minVal/maxVal for all dimensions.
if(length(numBin)<size(dataMat,3))
    numBin=repmat(numBin(1),[1,size(dataMat,3)]);
end
if(length(minVal)<size(dataMat,3))
    minVal=repmat(minVal(1),[1,size(dataMat,3)]);
end
if(length(maxVal)<size(dataMat,3))
    maxVal=repmat(maxVal(1),[1,size(dataMat,3)]);
end

% Compute bin indices
binInd=ones(size(dataMat,1),size(dataMat,2));
binOffset=[1 numBin];
for i=1:size(dataMat,3)
    tempInd=max(min(ceil(numBin(i)*(double(dataMat(:,:,i))-minVal(i))/(maxVal(i)-minVal(i))),numBin(i)),1);
    binInd=binInd+prod(binOffset(1:i))*(tempInd-1);
end

%% Matlab version
r=size(binInd,1); c=size(binInd,2);
intHist=zeros(r+1,c+1,prod(numBin),'int32');

intHist(2,2,binInd(1,1))=1;
for i=2:c
    intHist(2,i+1,:)=intHist(2,i,:);
    intHist(2,i+1,binInd(1,i))=intHist(2,i+1,binInd(1,i))+1;
end
for j=2:r
    intHist(j+1,2,:)=intHist(j,2,:);
    intHist(j+1,2,binInd(j,1))=intHist(j+1,2,binInd(j,1))+1;
end

for i=2:c
    for j=2:r
        intHist(j+1,i+1,:)=intHist(j,i+1,:)+intHist(j+1,i,:)-intHist(j,i,:);   
        intHist(j+1,i+1,binInd(j,i))=intHist(j+1,i+1,binInd(j,i))+1;
    end
end


