% Release
function [ v0reg,mreg,vx,vy,regionRows,regionCols, Err0,Err] = RegionBased( I1,I2,...
    RegionNum,GlobalCorrection,theta,DefRestrainTrans,method)
%参数：
% I1,I2: 输入图像
% RegionNum: 分块数,格式为[RowsNum ColsNum]
% GlobalCorrection: 全局修正,1为使用,0为不使用.默认为1.
% theta: Radon变换使用的角度.默认为1:180.
% DefRestrainTrans: 是否将估计模型简化为平移.1为是,0为否.默认0.
% method: 估计仿射变换时是否使用多分辨率. method.Global设置全局修正,method.Local设置分区估计.
%         默认为: method.Global=1; method.Local=0;
%
%返回值：
% v0reg: 1*2*N的三维数组,N为分区数,保存每个分区的v0分量.分区顺序为从左到右,从上到下.
% mreg: 类似v0reg,2*2*N的三维数组,保存每个分区的m分量.
% vx: 运动场x方向分量.
% vy: 运动场y方向分量.
% regionRows: 分区信息.行切割位置. [regionRows(i)+1 regionRows(i+1)]为第i行.
% regionCols: 分区信息.列切割位置. 同上.


    if nargin < 7
        method.Global=1;method.Local=0;
        if nargin < 6
            DefRestrainTrans = 0;
            if nargin < 5
                theta=1:180;
                if nargin < 4
                    GlobalCorrection = 1;
                    if nargin < 3
                        error(message('RegionBased:NotEnoughInputs'));
                    end
                end
            end
        end
    end
    
    if RegionNum(1)<=0 || RegionNum(2)<=0
        error('RegionNum error');
    end

    v0=[0 0];
    
    Err0=0;
    
    if GlobalCorrection==1        
        if method.Global==1
            [ v0,~,vxg,vyg,Err0 ] = Multiscale( I1,I2,5,theta,1);
        elseif method.Global==0
            [ v0,~,Err0] = EstimateAffine( I1,I2,theta,1);
            vxg = v0(1).*ones(size(I1));
            vyg = v0(2).*ones(size(I1));
        end
        
        A(1:2,1:2)=eye(2);
        A(3,1:2)=v0;

        picCenter=floor((size(I1)+1)/2);
        T=maketform('affine',A);
        I1pyd=imtransform(I1,T,'FillValues',0,'UData',[-picCenter(2)+1,size(I1,2)-picCenter(2)],'VData',[size(I1,1)-picCenter(1),-picCenter(1)+1],...
        'XData',[-picCenter(2)+1,size(I1,2)-picCenter(2)],'YData',[size(I1,1)-picCenter(1),-picCenter(1)+1],'Size',size(I1));

%         I1_ori=I1;
%         I2_ori=I2;
        if v0(1)<0
            I1=I1pyd(:,1:end+ceil(v0(1)));
            I2=I2(:,1:end+ceil(v0(1)));
        elseif v0(1)>0
            I1=I1pyd(:,1+ceil(v0(1)):end);
            I2=I2(:,1+ceil(v0(1)):end);
        end
    end
    
    regionRows=zeros(1,RegionNum(1)+1);regionCols=zeros(1,RegionNum(2)+1);
    regionRowsT=floor(size(I1,1)/RegionNum(1));regionColsT=floor(size(I1,2)/RegionNum(2));
    regionRows(1)=0;
    if RegionNum(1)==1
        regionRows(2)=size(I1,1);
    else
        for i=1:RegionNum(1)-1
            regionRows(i+1)=i*regionRowsT;
        end
        regionRows(i+2)=size(I1,1);
    end
    
    regionCols(1)=0;
    if RegionNum(1)==1
        regionCols(2)=size(I1,2);
    else
        for i=1:RegionNum(2)-1
            regionCols(i+1)=i*regionColsT;
        end
        regionCols(i+2)=size(I1,2);
    end

    v0reg=zeros(1,2,RegionNum(1)*RegionNum(2));
    mreg =zeros(2,2,RegionNum(1)*RegionNum(2));    
    vx=zeros(size(I1));vy=zeros(size(I1));
    Err=zeros(RegionNum(1),RegionNum(2));
    count=1;
    for i=1:RegionNum(1)
        for j=1:RegionNum(2)
            I1reg=I1(regionRows(i)+1:regionRows(i+1),regionCols(j)+1:regionCols(j+1));
            I2reg=I2(regionRows(i)+1:regionRows(i+1),regionCols(j)+1:regionCols(j+1));
            
            if method.Local==0
                [v0t,mt, Err(i,j)]=EstimateAffine(I1reg,I2reg,theta,DefRestrainTrans);
                
                picCenter=floor((size(I1reg)+1)/2);
                [X,Y]=meshgrid(-picCenter(2)+1:size(I1reg,2)-picCenter(2),...
                    size(I1reg,1)-picCenter(1):-1:-picCenter(1)+1);
                 vxreg=mt(1,1)*X+mt(1,2)*Y+v0t(1);
                 vyreg=mt(2,1)*X+mt(2,2)*Y+v0t(2);
                
            elseif method.Local==1
                [v0t,mt,vxreg,vyreg, Err(i,j)]=Multiscale( I1reg,I2reg,3,theta,DefRestrainTrans);
            end
            
            vx(regionRows(i)+1:regionRows(i+1),regionCols(j)+1:regionCols(j+1))=vxreg(:,:);
            vy(regionRows(i)+1:regionRows(i+1),regionCols(j)+1:regionCols(j+1))=vyreg(:,:);
            v0reg(:,:,count)=v0t+v0;
            mreg(:,:,count)=mt;
            
            count=count+1;
        end
    end
    
    if GlobalCorrection==1
%         I1=I1_ori;I2=I2_ori;
        if v0(1)<0
            vx=[zeros(size(I1,1),-ceil(v0(1))) vx];
            vy=[zeros(size(I1,1),-ceil(v0(1))) vy];
        elseif v0(1)>0
            vx=[vx zeros(size(I1,1),ceil(v0(1)))];
            vy=[vy zeros(size(I1,1),ceil(v0(1)))];
        end
        vx=vx+vxg;
        vy=vy+vyg;
    end
    
end
