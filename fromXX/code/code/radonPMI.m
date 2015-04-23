if ~exist('EVA','var')
    clear variables
end
close all
sigma=1;

OUTPUT=~exist('EVA','var');
DefZero=0.00001;
DefUseQ=1;
DefUseR=1;
DefMyRadon=0; %just for [0 45 90 135]
DefReduce=1;
DefExpandImg=0;
DefShowImgs=0;
DefRectNormlized=1;
DefIntensityNormlized=0;

DefShowPara=0;
Expand=20;
DefExhibit=OUTPUT&0;Print=1;
DefShowFinallDiff=0;

if ~exist('EVA','var')
DefNoisy=0;  Noise=0.03;%0.07 0.0018 0.005 0.018 Noise=10;
DefTransType=1; TransType=DefTransType;
DefExpand=1;
end

DefNoisyBefore=0;   
MultiHeight=3;
DefMode=1;
DefShowVelocity=0;DefEstErr=1;
DefShowGroundTruth=0;interval=50;
%theta = [0 atan(size(I1,1)/size(I1,2))/pi*180 90 180-atan(size(I1,1)/size(I1,2))/pi*180];%0:180;
%theta = [0 36 90 144];
theta = [0 45 90 135];
% theta = 1:30:180;

while(DefTransType==1)
    if ~exist('img1','var');
%     M=[0.05   0.01
%        0.01   0.06];
%     V=[0.5   0.5];
%   M=[0.05   0.010
%      0.080   0.04];
%   V=[10.0   20.0];
%    M=[-0.3   0.120
%     0.120   -0.5];
%   V=[10.0   20.0];
  M=[0.02   0.010
     0.20   0.04];
  V=[0.1   0.2];
%     M=[0.07   -0.12
%      0.09   0.1];
%   V=[1   2];
%    M=[0.1   0.120
%      0.090   -0.2];
%   V=[3   0.5];

%   M=[0.006 0.068
%     -0.034 0.012];
% V=[0.1 0.2];
%   img1='lena.jpg';
  %   img1='test3.jpg';
%     img1='forest3.bmp';
 %   img1='office.bmp';
end



    %I1=imread('new2binarytreed-03.tif');
    I1=rgb2gray(imread(img1));
%     I1=imresize(I1,[125,floor(size(I1,2)/size(I1,1)*125)]);

%     I1=zeros(161,201);
%     for i=-39:41
%         for j=-49:51
%             I1(80+i,100+j)=255;
%         end
%     end

%     I1=zeros(160,160);
%     for i=1:160
%             I1(i,80)=225;
%             I1(80,i)=255;
%     end


    if DefExhibit==1
        A(1:2,1:2)=eye(2)+M';
        A(3,1:2)=V;
%         I1o=flipud(I1);
        picCenter=floor((size(I1)+1)/2);
        T=maketform('affine',A);
        I2o=imtransform(I1,T,'FillValues',255,'UData',[-picCenter(2)+1,size(I1,2)-picCenter(2)],'VData',[size(I1,1)-picCenter(1),-picCenter(1)+1],...
    'XData',[-picCenter(2)+1,size(I1,2)-picCenter(2)],'YData',[size(I1,1)-picCenter(1),-picCenter(1)+1],'Size',size(I1));
%         I2o=flipud(I2o);
        I1o=I1;
        figure('Name','original image');
        imshow(I1o);
        axis off
        figure('Name','deformed image');
        imshow(I2o);
        axis off
    end

    if DefExpand==1
        I1tmp=zeros(size(I1)+2.*[Expand Expand]);
        I1tmp(Expand+1:end-Expand,Expand+1:end-Expand)=I1;
        I1=uint8(I1tmp);
    end
    
    w=fspecial('gaussian',[3 3],1);
      w=fspecial('average',[5 5]);
      I1=imfilter(I1,w);
      
    DefReduce=DefReduce&~DefMyRadon;
    DefExpandImg=DefExpandImg|DefMyRadon;

    !echo ==============================================
    fprintf('img:%s UseQ:%d UseR:%d Reduce:%d RectNorm:%d InteNorm:%d\nimgSize:[%d %d] TransType:1\n',...
            img1,DefUseQ,DefUseR,DefReduce,DefRectNormlized,DefIntensityNormlized,...
            size(I1,1),size(I1,2));

    flatI=ones(size(I1));
    
    if DefExpandImg==1
        if size(I1,1)<size(I1,2)
            addrows=size(I1,2)-size(I1,1);
            z1=zeros(ceil(addrows/2),size(I1,2));
            z2=zeros(addrows-ceil(addrows/2),size(I1,2));
            flatI=[z1;ones(size(I1));z2];
            Itmp=[z1;I1;z2];
            I1=Itmp;
        elseif size(I1,1)>size(I1,2)
            addcols=size(I1,1)-size(I1,2);
            z1=zeros(size(I1,1),ceil(addcols/2));
            z2=zeros(size(I1,1),addcols-ceil(addcols/2));
            flatI=[z1 ones(size(I1)) z2];
            Itmp=[z1 I1 z2];
            I1=Itmp;
        end
    end

%     A(1:2,1:2)=eye(2)-M;
%     A(1:2,3)=[0;0];
%     A(3,1:2)=[V(1,1) V(1,2)]; %debug!
%     A(3,3)=1;
% 
%     picCenter=floor((size(I1)+1)/2);
%     T=maketform('affine',A);
%     % I2tmp=imtransform(I1tmp,T,'UData',[-picCenter(2)+1,size(I1,2)-picCenter(2)],'VData',[-picCenter(1)+1,size(I1,1)-picCenter(1)],...
%     % 'XData',[-picCenter(2)+1,size(I1,2)-picCenter(2)],'YData',[-picCenter(1)+1,size(I1,1)-picCenter(1)],'FillValues',0);
%     I2=imtransform(I1,T,'FillValues',0,'UData',[-picCenter(2)+1,size(I1,2)-picCenter(2)],'VData',[-picCenter(1)+1,size(I1,1)-picCenter(1)],...
%     'XData',[-picCenter(2)+1,size(I1,2)-picCenter(2)],'YData',[-picCenter(1)+1,size(I1,1)-picCenter(1)]);
%     % I2=flipud(I2tmp);
%     %

    A(1:2,1:2)=eye(2)+M';
    A(3,1:2)=V;

    picCenter=floor((size(I1)+1)/2);
    T=maketform('affine',A);
    I2=imtransform(I1,T,'FillValues',0,'UData',[-picCenter(2)+1,size(I1,2)-picCenter(2)],'VData',[size(I1,1)-picCenter(1),-picCenter(1)+1],...
    'XData',[-picCenter(2)+1,size(I1,2)-picCenter(2)],'YData',[size(I1,1)-picCenter(1),-picCenter(1)+1],'Size',size(I1));


    if DefNoisy==1
        %     I2_ori=I2;
        I1tmp=I1;I2tmp=I2;
        I1=imnoise(I1,'gaussian',0,Noise);
        I2=imnoise(I2,'gaussian',0,Noise);

        % I1=awgn(double(I1),Noise,'measured');I2=awgn(double(I2),Noise,'measured');

        % I1=AddGaussianNoise(I1,20);I2=AddGaussianNoise(I2,20);
        snr1=SNR(double(I1tmp),double(I1));snr2=SNR(double(I2tmp),double(I2));
        if Print~=0
            fprintf('SNR1:%f SNR2:%f\n',snr1,snr2);
        end
    end
    
%     figure(12);
%     imshow(I1);
%     H=fspecial('gaussian');
%     I1=imfilter(I1,H);
%     I2=imfilter(I2,H);

    DefTransType=0;
end
while(DefTransType==2)
    img1='new2binarytreed-20.tif';
    img2='new2binarytreed-21.tif';
%     
%     img1='left.ppm';
%     img2='right.ppm';
    
    I1=(imread(img1));
    I2=(imread(img2));
%     I1=imresize(I1,[size(I1,1)/4,size(I1,2)/4]);
%     I2=imresize(I2,[size(I2,1)/4,size(I2,2)/4]);
    flatI=ones(size(I1));

    DefReduce=DefReduce&~DefMyRadon;
    DefExpandImg=DefExpandImg|DefMyRadon;

    !echo ==============================================
    fprintf('img1:%s img2:%s\nUseQ:%d UseR:%d Reduce:%d RectNorm:%d InteNorm:%d\nimgSize:[%d %d] TransType:2\n',...
            img1,img2,DefUseQ,DefUseR,DefReduce,DefRectNormlized,DefIntensityNormlized,...
            size(I1,1),size(I1,2));

    if DefExpandImg==1
        if size(I1,1)<size(I1,2)
            addrows=size(I1,2)-size(I1,1);
            z1=zeros(ceil(addrows/2),size(I1,2));
            z2=zeros(addrows-ceil(addrows/2),size(I1,2));
            flatI=[z1;ones(size(I1));z2];
            Itmp=[z1;I1;z2];
            I1=Itmp;
        elseif size(I1,1)>size(I1,2)
            addcols=size(I1,1)-size(I1,2);
            z1=zeros(size(I1,1),ceil(addcols/2));
            z2=zeros(size(I1,1),addcols-ceil(addcols/2));
            flatI=[z1 ones(size(I1)) z2];
            Itmp=[z1 I1 z2];
            I1=Itmp;
        end
    end
            I1tmp=I1;I2tmp=I2;
        I1=imnoise(I1,'gaussian',0,Noise);
        I2=imnoise(I2,'gaussian',0,Noise);

        % I1=awgn(double(I1),Noise,'measured');I2=awgn(double(I2),Noise,'measured');

        % I1=AddGaussianNoise(I1,20);I2=AddGaussianNoise(I2,20);
        snr1=SNR(double(I1tmp),double(I1));snr2=SNR(double(I2tmp),double(I2));
            fprintf('SNR1:%f SNR2:%f\n',snr1,snr2);
        
    DefTransType=0;
end

height=size(I1,1);
width=size(I1,2);

if DefIntensityNormlized==1
totalI1=sum(sum(I1));
totalI2=sum(sum(I2));
I2=I2./(totalI2/totalI1);
end

if DefShowImgs==1
    figure(10);
    imshow(I1);
    figure(11);
    imshow(I2);
%     figure(13);
%     imshow(abs(I2-I1));
% figure(3);
% imshow(I1tmp);
% figure(4);
% imshow(I2tmp);
end

% thetatest=0:180;
% [Radon1,xp1]=radon(I1,thetatest);
% [Radon2,xp2]=radon(I2,thetatest);
% figure(25);
% imshow(Radon1,[],'Xdata',thetatest,'Ydata',xp1,'InitialMagnification','fit');
% xlabel('\theta (degrees)')
% ylabel('x''')
% colormap(hot), colorbar
% iptsetpref('ImshowAxesVisible','off')
% figure(26);
% imshow(Radon2,[],'Xdata',thetatest,'Ydata',xp1,'InitialMagnification','fit');
% xlabel('\theta (degrees)')
% ylabel('x''')
% colormap(hot), colorbar
% iptsetpref('ImshowAxesVisible','off')
% figure(27);
% imshow(Radon2-Radon1,[],'Xdata',thetatest,'Ydata',xp1,'InitialMagnification','fit');
% xlabel('\theta (degrees)')
% ylabel('x''')
% colormap(hot), colorbar
% iptsetpref('ImshowAxesVisible','off')
tic
if DefMode==0
    [v0 m]=EAPMI(I1,I2,theta,img1);
else
    [ v0,m,vx,vy ] = MultiscalePMI(I1,I2,MultiHeight-1,theta,img1,0);
end
toc
disp(' ');
disp('v0 = ');
disp(v0);
disp('m = ');
disp(m);

if DefShowFinallDiff==1 && (DefMode==0||DefMode==1) 
    if DefExhibit==1
        A(1:2,1:2)=(eye(2)+m')^-1;
        A(3,1:2)=-v0*A(1:2,1:2);
        picCenter=floor((size(I1o)+1)/2);
        T=maketform('affine',A);
        I3=imtransform(I2o,T,'FillValues',0,'UData',[-picCenter(2)+1,size(I1o,2)-picCenter(2)],'VData',[size(I1o,1)-picCenter(1),-picCenter(1)+1],...
        'XData',[-picCenter(2)+1,size(I1o,2)-picCenter(2)],'YData',[size(I1o,1)-picCenter(1),-picCenter(1)+1],'Size',size(I1o));
    else
        A(1:2,1:2)=(eye(2)+m')^-1;
        A(3,1:2)=-v0*A(1:2,1:2);
        picCenter=floor((size(I1)+1)/2);
        T=maketform('affine',A);
        I3=imtransform(I2,T,'FillValues',0,'UData',[-picCenter(2)+1,size(I1,2)-picCenter(2)],'VData',[size(I1,1)-picCenter(1),-picCenter(1)+1],...
        'XData',[-picCenter(2)+1,size(I1,2)-picCenter(2)],'YData',[size(I1,1)-picCenter(1),-picCenter(1)+1],'Size',size(I1));
    end
    figure;
    imshow(abs(I1o-I3),[]);
%     title('diff');
end

%===============================================================================
    picCenter=floor((size(I1)+1)/2);
    [X,Y]=meshgrid(-picCenter(2)+1:size(I1,2)-picCenter(2),...size(I1,2)/10
            size(I1,1)-picCenter(1):-1:-picCenter(1)+1);
    vxc=zeros(size(X));vyc=zeros(size(X));
    
    if  DefMode==0
        vx=zeros(size(X));vy=zeros(size(X));
        for i=1:size(X,1)
            for j=1:size(X,2)
                v=v0'+m*[X(i,j) Y(i,j)]';
                vx(i,j)=v(1);  vy(i,j)=v(2);
            end
        end
    end
    
    if TransType==1        
        for i=1:size(X,1)
            for j=1:size(X,2)
                v=V'+M*[X(i,j) Y(i,j)]';
                vxc(i,j)=v(1); vyc(i,j)=v(2);
            end
        end

    elseif TransType==2 && strncmp(img1,'new2binarytreed',14)      
        ff=fopen('Div20_ux_150x150_float32.raw');
        raw=fread(ff,'*float32');
%         vxc=zeros(150,150);
        for i=1:150
            for j=1:150
                vxc(i,j)=raw((i-1)*150+j);
            end
        end
        ff=fopen('Div20_uy_150x150_float32.raw');
        raw=fread(ff,'*float32');
%         vyc=zeros(150,150);
        for i=1:150
            for j=1:150
                vyc(i,j)=raw((i-1)*150+j);
            end
        end
        
   elseif TransType==2 && strncmp(img1,'left',4)
        s=load('final_labels.mat');
        final_labels=s.final_labels;
%         vxc=zeros(size(I1));
        for i=1:size(I1,1)
            for j=1:size(I1,2)
                vxc(i,j)=-final_labels(i,j);
            end
        end
%         vyc=zeros(size(I1));
%         for i=1:size(I1,1)
%             for j=1:size(I1,2)
%                 vyc(i,j)=0;
%             end
%         end
    elseif TransType==2 && strncmp(img1,'yosemite',8)
        s=load('yosemite.mat');
        vxc=s.flow.vx;
        vyc=s.flow.vy;
    end
%===============================================================================

if DefShowVelocity==1

    X_disp=X(1:interval:end,1:interval:end);
    Y_disp=Y(1:interval:end,1:interval:end);

    vx_disp=vx(1:interval:end,1:interval:end);
    vy_disp=vy(1:interval:end,1:interval:end);
    
    iptsetpref('ImshowAxesVisible','on')
    figure;
    
    if DefExpand==1 && DefExhibit==1
%         I1show=I1(Expand+1:end-Expand,Expand+1:end-Expand);
        picCenterShow=floor((size(I1o)+1)/2);
        imshow(I1o,[],'XData',[-picCenterShow(2)+1,size(I1o,2)-picCenterShow(2)],'YData',[size(I1o,1)-picCenterShow(1),-picCenterShow(1)+1]);
    else
        imshow(I1,[],'XData',[-picCenter(2)+1,size(I1,2)-picCenter(2)],'YData',[size(I1,1)-picCenter(1),-picCenter(1)+1]);
    end
    axis xy
    hold on
    quiver(X_disp,Y_disp,vx_disp,vy_disp,'r');
   
    if DefMode==2
        for i=2:size(regionRows,2)-1
            plot(-picCenter(2)+1:size(I1,2)-picCenter(2),regionRows(i)-picCenter(1)+1,'-g');
        end
        for i=2:size(regionCols,2)-1
             plot(regionCols(i)-picCenter(2)+1,-picCenter(2)+1:size(I1,2)-picCenter(2),'-g');
        end
    end    
    
    if DefShowGroundTruth==1
        vxc_disp=vxc(1:interval:end,1:interval:end);
        vyc_disp=vyc(1:interval:end,1:interval:end);
        quiver(X_disp,Y_disp,vxc_disp,vyc_disp,'-g');
    end
    
    hold off

end

if DefEstErr==1

    ang=0;mag=0;
    for i=1:size(X,1)
        for j=1:size(X,2)
            Ve=[vx(i,j) vy(i,j) 1];
            Vc=[vxc(i,j) vyc(i,j) 1];
            ang=ang+acos(Ve*Vc'/(norm(Vc,2)*norm(Ve,2)));
            mag=mag+abs(norm(Vc,2)-norm(Ve,2));
        end
    end        
    ang=180*ang/(size(X,1)*size(X,2))/pi;
    mag=mag/(size(X,1)*size(X,2));
    
    fprintf('ang=%f  mag=%f\n',ang,mag);
end

