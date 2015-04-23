% Release
function [ v0, m ,Err ] = EstimateAffine( I1 , I2, theta, DefRestrainTrans)
% tic
%ESTIMATEAFFINE Summary of this function goes here
%   Detailed explanation goes here
    DefZero=0.00001;
    
%     if size(I1, 3)~=1
%          I1 = rgb2gray(I1);
%     end
%     if size(I2, 3)~=1
%          I2 = rgb2gray(I2);
%     end
%     
    if size(I1)~=size(I2)
        error('size of the two images should be the same');
    end
    
    flatI=ones(size(I1));    
    picCenter=floor((size(I1)+1)/2);
    
    [Radon1,xp1]=radon(I1,theta);
    [Radon2,~]=radon(I2,theta);
    [Radonflat,~]=radon(flatI,theta);

    Radon1=Radon1./Radonflat;
    Radon2=Radon2./Radonflat;
        Radon1(isnan(Radon1)==1)=0;
        Radon2(isnan(Radon2)==1)=0;
    
    theta1=cos(theta*pi/180);
    theta2=sin(theta*pi/180);
    coreGy=[0.5;0;-0.5];
    gp=conv2(coreGy,Radon1);
    gp=gp(2:size(gp,1)-1,:);

    %derivative of g with respect to theta====================================
    % x*f(x,y), y*f(x,y)----------------------------------------------
    Itheta1=zeros(size(I1));Itheta2=zeros(size(I1));
    for jj=1:size(I1,1)
        Itheta1(jj,:)=(-picCenter(2)+1:size(I1,2)-picCenter(2)).*double(I1(jj,:));
    end
    for jj=1:size(I1,2)
        Itheta2(:,jj)=-1.*(-picCenter(1)+1:size(I1,1)-picCenter(1))'.*double(I1(:,jj));
    end
%     D=diag(-picCenter(1)+1:size(I1,1)-picCenter(1));
%     Itheta2=-D*double(I1);
%     D=diag(-picCenter(2)+1:size(I1,2)-picCenter(2));
%     Itheta1=double(I1)*D;
    
    % radon of x*f(x,y), y*f(x,y)-----------------------------------
    RItheta1=radon(Itheta1,theta);
    RItheta2=radon(Itheta2,theta);

    % normalized
    RItheta1=RItheta1./Radonflat;
    RItheta2=RItheta2./Radonflat;
        RItheta1(isnan(RItheta1)==1)=0;
        RItheta2(isnan(RItheta2)==1)=0;

    gtheta1=-conv2(coreGy,RItheta1);
    gtheta2=-conv2(coreGy,RItheta2);
    gtheta1=gtheta1(2:size(gtheta1,1)-1,:);
    gtheta2=gtheta2(2:size(gtheta2,1)-1,:);

    Utheta=theta;
    Up=xp1;
    cols=size(Up,1)*size(Utheta,2);
    H=zeros(cols,6);R=zeros(cols,1);
    Q=zeros(cols);%spdiags
    count=1;
    for i=1:size(Utheta,2)
        xi=i;

         counti=1;
            while abs(Radonflat(counti,xi)-0)<DefZero
                counti=counti+1;
            end
            start=counti+10;

            counti=size(Radonflat,1);
            while abs(Radonflat(counti,xi)-0)<DefZero
                counti=counti-1;
            end
            ending=counti-10;
%         start=find(Radonflat(:,xi),1,'first')+10;
%         ending=find(Radonflat(:,xi),1,'last')-10;

        for j=start:ending;        
            yi=find(xp1==Up(j),1);

            a1=theta1(xi)*theta2(xi)^2*gtheta1(yi,xi)+Radon1(yi,xi)*theta2(xi)^2-...
                xp1(yi)*theta1(xi)^2*gp(yi,xi)-theta1(xi)^2*theta2(xi)*gtheta2(yi,xi);
            b1=theta1(xi)^3*gtheta2(yi,xi)-theta1(xi)^2*theta2(xi)*gtheta1(yi,xi)-...
                theta1(xi)*theta2(xi)*xp1(yi)*gp(yi,xi)-Radon1(yi,xi)*theta1(xi)*theta2(xi);
            c1=theta2(xi)^3*gtheta1(yi,xi)-theta2(xi)^2*theta1(xi)*gtheta2(yi,xi)-...
                theta1(xi)*theta2(xi)*xp1(yi)*gp(yi,xi)-Radon1(yi,xi)*theta1(xi)*theta2(xi);
            d1=theta2(xi)*theta1(xi)^2*gtheta2(yi,xi)+Radon1(yi,xi)*theta1(xi)^2-...
                xp1(yi)*theta2(xi)^2*gp(yi,xi)-theta2(xi)^2*theta1(xi)*gtheta1(yi,xi);
            e1=-theta1(xi)*gp(yi,xi);
            f1=-theta2(xi)*gp(yi,xi);
            
            H(count,:)=[a1 b1 c1 d1 e1 f1];
            R(count,1)=(Radon2(yi,xi)-Radon1(yi,xi))*1;
            
%             if(abs(Radonflat(yi,xi)-0)<DefZero)
%                 Q(count,count)=1;
%             else
                Q(count,count)=1/Radonflat(yi,xi);
%             end
        
            count=count+1;
        end   
    end

    if DefRestrainTrans==1
        H=H(1:count-1,5:6);
    else
        H=H(1:count-1,:);
    end
    
    R=R(1:count-1,:);
    Q=Q(1:count-1,1:count-1);
    res=H\R;
    if count~=1
        res=(H'/Q*H)\H'/Q*R;
    else
        res=H\R;
    end
%     res=H\R;
    Err=norm(H*res-R,2)+0.5*res(1);
    if DefRestrainTrans==1
        v0=res';
        m=[0 0;0 0];
    else
        v0=res(5:6,1)';
        m=[res(1:2,1)';res((3:4),1)'];
    end
% toc
end

