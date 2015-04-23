function [intHist,binInd]=formLBPIntegralHistogram(dataMat)
% currently only 1 histograms

% Compute bin indices
%binInd=ceil(numBin*(double(dataMat)+1e-5)/maxVal);
%binInd(binInd>numBin)=numBin;

if isempty(dataMat)
    error('Datamatrix is empty');
end

binInd = dataMat+1;

%% Matlab version
[r, c] = size(dataMat);

intHist=zeros(r+1,c+1,59,'int32');

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