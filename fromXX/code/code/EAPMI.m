function [ v0 m ] = EAPMI( I1,I2,theta,img1 )
%EAPMI Summary of this function goes here
%   Detailed explanation goes here
DefMyRadon=0;DefRectNormlized=1; DefReduce=1;DefUseQ=1;DefUseR=1;DefShowPara=0;
sigma=1;
DefZero=0.00001;
Psi=zeros(2,size(theta,2));
W=zeros(size(theta,2),2);
A=zeros(size(theta,2),3);
R=zeros(size(theta,2)*2,size(theta,2)*2);
flatI=ones(size(I1));
for i=1:size(theta,2)
    clear R1 xp1 xp2 Rflat xpflat R2 dR1;
    clear g1 g2 gp p H1 gt;
    tmpDefMyRadon=DefMyRadon;
    
    while(tmpDefMyRadon==1)
        if size(theta,2)~=4||theta(1)~=0||theta(2)~=45||theta(3)~=90||theta(4)~=135
            error('illegal thetas');
        end
        if theta(i)~=90
            plength=ceil(2*(width/2+height/2*abs(tan(theta(i)*pi/180))));%ceil(2*abs(cos(theta(i)*pi/180-atan(height/width)))*sqrt(height^2+width^2)/2/sqrt(2));
        else 
            plength=height;
        end
        R1=zeros(plength,1);xp1=zeros(plength,1);
        R2=zeros(plength,1);xp2=zeros(plength,1);
        Rflat=zeros(plength,1);xpflat=zeros(plength,1);
        if theta(i)==0
            for jj=1:width
                R1(jj)=sum(I1(:,jj));xp1(jj)=jj-1-picCenter(2);
                R2(jj)=sum(I2(:,jj));xp2(jj)=jj-1-picCenter(2);
                Rflat(jj)=height;xpflat(jj)=jj-1-picCenter(2);
            end            
        elseif theta(i)==45
                for jj=1:plength
                    sum1=0;sum2=0;sum3=0;
                    if jj<=height
                        for ii=0:jj-1
                            sum1=sum1+I1(ii+height-jj+1,ii+1);
                            sum2=sum2+I2(ii+height-jj+1,ii+1);
                            sum3=sum3+flatI(ii+height-jj+1,ii+1);
                        end
                    else
                        for ii=0:width+height-jj
                            sum1=sum1+I1(ii+1,jj-height+ii);
                            sum2=sum2+I2(ii+1,jj-height+ii);
                            sum3=sum3+flatI(ii+1,jj-height+ii);                            
                        end
                    end
                    R1(jj)=sum1;xp1(jj)=(jj-height)/sqrt(2);
                    R2(jj)=sum2;xp2(jj)=(jj-height)/sqrt(2);%(jj-1-picCenter(2))/sqrt(2)
                    Rflat(jj)=sum3;xpflat(jj)=(jj-height)/sqrt(2);
                end
        elseif theta(i)==90
            for jj=1:height
                R1(jj)=sum(I1(jj,:));xp1(jj)=jj-1-picCenter(1);
                R2(jj)=sum(I2(jj,:));xp2(jj)=jj-1-picCenter(1);
                Rflat(jj)=width;xpflat(jj)=jj-1-picCenter(1);
            end
        elseif theta(i)==135
            for jj=1:plength
                sum1=0;sum2=0;sum3=0;
                    if jj<=height
                        for ii=0:jj-1
                            sum1=sum1+I1(jj-ii,ii+1);
                            sum2=sum2+I2(jj-ii,ii+1);
                            sum3=sum3+flatI(jj-ii,ii+1);
                        end
                    else
                        for ii=0:2*height-jj
                            sum1=sum1+I1(jj-height+ii,width-ii);
                            sum2=sum2+I2(jj-height+ii,width-ii);
                            sum3=sum3+flatI(jj-height+ii,width-ii);                            
                        end
                    end
                    R1(jj)=sum1;xp1(jj)=(jj-height)/sqrt(2);%(jj-1-picCenter(2))/sqrt(2);
                    R2(jj)=sum2;xp2(jj)=(jj-height)/sqrt(2);
                    Rflat(jj)=sum3;xpflat(jj)=(jj-height)/sqrt(2);
            end
        end
        tmpDefMyRadon=-1;
    end
    
    while(tmpDefMyRadon==0)
        [R1 xp1]=radon(I1,theta(i));
        [R2 ~]=radon(I2,theta(i));
        [Rflat ~]=radon(flatI,theta(i));
        
        tmpDefMyRadon=-1;
    end

%     %for debug
%     flatInPaper=zeros(size(xp1,1),1);

if DefRectNormlized==1
    R1=R1./Rflat;
    R2=R2./Rflat;
    R1(isnan(R1)==1)=0;
    R2(isnan(R2)==1)=0;
end
    
    %for debug
%     Radon1=[R1 xp1];
%     Radon2=[R2 xp2];
%     Radonflat=[Rflat xpflat];

    covR1=conv(R1,[1 0 -1]');
    covxp1=conv(xp1,[1 0 -1]');
    dR1=[0;(covR1(3:(size(covR1,1)-2)))./(covxp1(3:(size(covxp1,1)-2)));0];
    %dxp1=xp1(2:(size(xp1,1)-1));

    %dR1 =[diff(R1)./diff(xp1);0];
    %dR1 =[0;diff(R1)./diff(xp1)];
    
    if DefReduce==1
        
        counti=1;
        while abs(Rflat(counti)-0)<DefZero
            counti=counti+1;
        end
%         if counti==1
%             start=2;
%         else
            start=counti+1;
%         end
        
        counti=size(Rflat,1);
        while abs(Rflat(counti)-0)<DefZero
            counti=counti-1;
        end
%         if counti==size(Rflat,1);
%             ending=size(Rflat,1)-1;
%         else
            ending=counti-1;
%         end

       % start=13; ending=202;
        g1=R1(start:ending);
        g2=R2(start:ending);
        gf=Rflat(start:ending);
        gp=dR1(start:ending);
        p=xp1(start:ending);
    else
        g1=R1;
        g2=R2;
        gf=Rflat;
        %dR1(abs(dR1)<0.0000000001)=0;
        gp=dR1;%[0;dR1];
        p=xp1;
    end
    
    %for debug
%     gpp=[gp p];
    
    H1=[gp gp.*p];
    gt=(g2-g1);    

%QQQQQQQQQQQQQQQQQQQQQQQQQ
    Q=zeros(size(p,1));
    
%     theta1=abs(sin(theta(i)*pi/180));
%     theta2=abs(cos(theta(i)*pi/180));    
%     for j=1:size(p,1)
%         if (abs(theta1-0)<DefZero)
% 				t=height;
%         elseif(abs(theta2-0)<DefZero)
% 				t=width;
%         else
% 				t=min(p(j)*theta2/theta1+width/2/theta1,-p(j)*theta1/theta2+height/2/theta2)-...
%                     max(p(j)*theta2/theta1-width/2/theta1,-p(j)*theta1/theta2-height/2/theta2);
%         end
%         Q(j,j)=sigma^2/(t);
        
%         %for debug
%         flatInPaper(j)=t;
%     end

    for j=1:size(p,1)
        if(abs(gf(j)-0)<DefZero)
            Q(j,j)=1;
        else
            Q(j,j)=sigma^2/gf(j);
        end
    end
%QQQQQQQQQQQQQQQQQQQQQQQQQ

    %for debug
%     eqt=[gp gp.*p gt];
    
    if DefUseQ==1
        Psi(:,i)=-(H1'/Q*H1)\H1'/Q*gt;
    else
        Psi(:,i)=-H1\gt;
    end
    %Psi(:,i)=-(H1'*H1)\H1'*gt;
    %Psi(:,i)=-inv((H1'*inv(Q)*H1))*H1'*inv(Q)*gt;
    
    CovPsi=inv(H1'/Q*H1);
    R(i,i)=CovPsi(1,1);
    R(i,i+size(theta,2))=CovPsi(1,2);
    R(i+size(theta,2),i)=CovPsi(2,1);
    R(i+size(theta,2),i+size(theta,2))=CovPsi(2,2);
    
    W(i,:)=[cos(theta(i)*pi/180) sin(theta(i)*pi/180)];
    A(i,:)=[cos(theta(i)*pi/180)^2 sin(theta(i)*pi/180)^2 2*cos(theta(i)*pi/180)*sin(theta(i)*pi/180)];
    
    if DefShowPara==1
    figure(i);
    plot(p,gp,p,gp.*p/100,p,gt,p,-Psi(1,i)*gp-Psi(2,i)*gp.*p,'DisplayName','H1');
    title([num2str(theta(i)) ' - ' img1]);
    legend('gp','gpp/100','gt','gtcal');
    axis tight ;
    xlabel('p');
    end
end 

if DefUseR==1
    D=[W zeros(size(W,1),size(A,2));zeros(size(A,1),size(W,2)) A];
    y=[Psi(1,:)';Psi(2,:)'];
    phi=(D'/R*D)\D'/R*y;
    v0=[phi(1) phi(2)];%debug!
    m =[phi(3) phi(5)/2;phi(5)/2 phi(4)]; %debug!
else
    v0= (W\Psi(1,:)')';
    mtmp = A\Psi(2,:)';
    m =[mtmp(1) mtmp(3)/2;mtmp(3)/2 mtmp(2)];%debug!
    % v0= (Wi'*Wi)^-1*Wi'*alphau0(1,:)'
    % m = (Ai'*Ai)^-1*Ai'*alphau0(2,:)'
end

end

