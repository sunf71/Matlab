clear all;
i1 = imread('in000007.jpg');
i2 = imread('in000008.jpg');
I1 = rgb2gray(i1);
I2 = rgb2gray(i2);
imwrite(I1,'i1.jpg');
imwrite(I2,'i2.jpg');
%# Create the gaussian filter with hsize = [5 5] and sigma = 2
G = fspecial('gaussian',[5 5],2);
%# Filter it
I1f = imfilter(I1,G,'same');
I2f = imfilter(I2,G,'same');
[v0,m,vx,vy]=Multiscale(I1f,I2f,10,0:10:179,0);
 A(1:2,1:2)=(eye(2)+m');
 A(3,1:2)=v0*A(1:2,1:2);
 picCenter=floor((size(i1)+1)/2);
 T=maketform('affine',A);
 I3=imtransform(I1,T,'FillValues',0,'UData',[-picCenter(2)+1,size(I1,2)-picCenter(2)],'VData',[size(I1,1)-picCenter(1),-picCenter(1)+1],...
        'XData',[-picCenter(2)+1,size(I1,2)-picCenter(2)],'YData',[size(I1,1)-picCenter(1),-picCenter(1)+1],'Size',size(I1));
 figure;
 imshow(I3,[]);
 imwrite(I3,'i3.jpg');
 figure;
 Id = abs(I2 - I3);
 imshow(Id,[]);
 imwrite(Id,'i2-i3.jpg');
