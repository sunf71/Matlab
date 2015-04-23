% Release
function [ v0,m,vx,vy, Err ] = Multiscale(I1,I2,height,theta,DefRestrainTrans)

    if height<0
        error('Multiscale height error');
    end
    mh=[0 0;0 0];
    v0h=[0 0];
    pyd=2;
    while height>=0

        I1pyd=imresize(I1,size(I1)/pyd^height);
        I2pyd=imresize(I2,size(I2)/pyd^height);

        A(1:2,1:2)=(eye(2)+mh');
        A(3,1:2)=pyd*[v0h(1,1) v0h(1,2)];

        picCenter=floor((size(I1pyd)+1)/2);
        T=maketform('affine',A);
        I1pyd=imtransform(I1pyd,T,'FillValues',0,'UData',[-picCenter(2)+1,size(I1pyd,2)-picCenter(2)],'VData',[size(I1pyd,1)-picCenter(1),-picCenter(1)+1],...
        'XData',[-picCenter(2)+1,size(I1pyd,2)-picCenter(2)],'YData',[size(I1pyd,1)-picCenter(1),-picCenter(1)+1],'Size',size(I1pyd));

        [v0r,mr, Err]=EstimateAffine(I1pyd,I2pyd,theta,DefRestrainTrans);
        
        mh=(eye(2)+mh)*(eye(2)+mr)-eye(2);
%         v0h=pyd*v0h*(eye(2)+mr)+v0r;
        v0h=pyd*v0h+v0r;
        
        height=height-1;
    end

    v0=v0h;
    m=mh;
        
    picCenter=floor((size(I1)+1)/2);
    [X,Y]=meshgrid(-picCenter(2)+1:size(I1,2)-picCenter(2),...
        size(I1,1)-picCenter(1):-1:-picCenter(1)+1);

        vx=m(1,1)*X+m(1,2)*Y+v0(1);
        vy=m(2,1)*X+m(2,2)*Y+v0(2);
end
