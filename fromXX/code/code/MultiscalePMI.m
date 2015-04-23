function [ v0,m,vx,vy ] = MultiscalePMI(I1,I2,height,theta,img1,MultiDebug)

    if height<0
        error('Multiscale height error');
    end
    ifigure=1;
    mh=[0 0;0 0];
    v0h=[0 0];
    isfirst=1;
    pyd=2;
    while height>=0
%         hgausspymd = vision.Pyramid;
%         hgausspymd.PyramidLevel = height;
%         I1pyd=mat2gray(double(step(hgausspymd, I1)));
%         I2pyd=mat2gray(double(step(hgausspymd, I2)));
%         release(hgausspymd);

        I1pyd=imresize(I1,size(I1)/pyd^height);
        I2pyd=imresize(I2,size(I2)/pyd^height);
%         w=fspecial('gaussian',[5 5],3);
%         I1pyd=imfilter(I1pyd,w);
%         I2pyd=imfilter(I2pyd,w);
    
        if MultiDebug==1
            figure(ifigure);
            ifigure=ifigure+1;
            imshow(I1pyd);
            title(['h=' num2str(height) '- I1']);
            figure(ifigure);
            ifigure=ifigure+1;
            imshow(I2pyd);
            title(['h=' num2str(height) '- I2 ori']);
            figure(ifigure);
            ifigure=ifigure+1;
            imshow(abs(I1pyd-I2pyd),[0,255]);
            title(['h=' num2str(height) '- diff before warped']);
        end
        
        A(1:2,1:2)=(eye(2)+mh');
        A(3,1:2)=pyd*[v0h(1,1) v0h(1,2)];

        picCenter=floor((size(I1pyd)+1)/2);
        T=maketform('affine',A);
        I1pyd=imtransform(I1pyd,T,'FillValues',0,'UData',[-picCenter(2)+1,size(I1pyd,2)-picCenter(2)],'VData',[size(I1pyd,1)-picCenter(1),-picCenter(1)+1],...
        'XData',[-picCenter(2)+1,size(I1pyd,2)-picCenter(2)],'YData',[size(I1pyd,1)-picCenter(1),-picCenter(1)+1],'Size',size(I1pyd));
        
%     cut=floor(size(I1pyd)*0.05/2);
%     I1pyd=I1pyd(cut(1):size(I1pyd,2)-cut(1),cut(2):size(I1pyd,2)-cut(2));
%     I2pyd=I2pyd(cut(1):size(I2pyd,2)-cut(1),cut(2):size(I2pyd,2)-cut(2));
%     
        if MultiDebug==1
            if isfirst==1
                isfirst=0;
            else
                figure(ifigure);
                ifigure=ifigure+1;
                imshow(I1pyd);
                title(['h=' num2str(height) '- I1 warped by h=' num2str(height+1) ' estimate']);

                figure(ifigure);
                ifigure=ifigure+1;
                imshow(abs(I1pyd-I2pyd),[0,255]);
                title(['h=' num2str(height) '- diff after warped by h=' num2str(height+1) ' estimate']);
            end
        end
        
        [v0r,mr]=EAPMI(I1pyd,I2pyd,theta,img1);
        
%         mh=mh+mr;
        mh=(eye(2)+mh)*(eye(2)+mr)-eye(2);
%         v0h=pyd*v0h*(eye(2)+mr)+v0r;
        v0h=pyd*v0h+v0r;
        
        if MultiDebug==1
            fprintf('h=%d \nv0r=[%g %g] mr=[%g %g %g %g]\nv0h=[%g %g] mh=[%g %g %g %g]\n',...
                height,v0r(1),v0r(2),mr(1,1),mr(1,2),mr(2,1),mr(2,2),...
                v0h(1),v0h(2),mh(1,1),mh(1,2),mh(2,1),mh(2,2));
        end
        
        height=height-1;
    end
    
    if MultiDebug==1
        A(1:2,1:2)=(eye(2)+mh');
        A(3,1:2)=[v0h(1,1) v0h(1,2)];

        picCenter=floor((size(I2)+1)/2);
        T=maketform('affine',A);
        I1pyd=imtransform(I1,T,'FillValues',0,'UData',[-picCenter(2)+1,size(I1,2)-picCenter(2)],'VData',[size(I1,1)-picCenter(1),-picCenter(1)+1],...
        'XData',[-picCenter(2)+1,size(I1,2)-picCenter(2)],'YData',[size(I1,1)-picCenter(1),-picCenter(1)+1],'Size',size(I1));

        I2pyd=I2;

        figure(ifigure);
        ifigure=ifigure+1;
        imshow(uint8(I1pyd));
        title(['finall - I1 warped by h=' num2str(height+1) ' estimate']);
        figure(ifigure);
        imshow(abs(I1pyd-I2pyd),[0,255]);
        title(['finall - diff after warped by h=' num2str(height+1) ' estimate']);
    end
    
    v0=v0h;
    m=mh;
    
%     if debug.DefShowImgs==1 
%         iptsetpref('ImshowAxesVisible','off');
%         figure(10);
%         imshow(I1,[]);%,'XData',[160 1],'YData',[160 1]        
%         title('I1');
%         axis image
%         figure(11);    
%         imshow(I2,[]); 
%         title('I2');
%         axis image
%     end
    
    picCenter=floor((size(I1)+1)/2);
    [X,Y]=meshgrid(-picCenter(2)+1:size(I1,2)-picCenter(2),...size(I1,2)/10
        size(I1,1)-picCenter(1):-1:-picCenter(1)+1);%size(I1,1)/10
    vx=zeros(size(X));vy=zeros(size(X));
    for i=1:size(X,1)
        for j=1:size(X,2)
            v=[v0(1);v0(2)]+m*[X(i,j) Y(i,j)]';
            vx(i,j)=v(1);
            vy(i,j)=v(2);
        end
    end
end

