文件说明：
    主要文件为Evaluation.m、main.m、SiftBased.m以及radonPMI.m。其中Evaluation.m为批量测试入口，论文中大多数实验数据通过这个文件得到；main.m是我们的算法实现的主文件，其中调用的Multiscale.m实现多分辨率计算，RegionBased.m实现分块计算，EstimateAffine.m实现核心算法；SiftBased.m是基于SIFT的配准算法的主文件，其中调用了siftWin32.exe计算SIFT特征，以及Ransc目录下的一些辅助函数；radonPMI.m为[20]中算法实现的主文件，其中调用了MultiscalePMI.m实现多分辨率计算，以及EAPMI.m实现核心算法。另有一些如图像、GroundTrueth数据等资源文件。

实验参数：
1.Diverging Tree实验
    Evaluation.m中，使用%% nosie部分，注释%% 8 images部分，并设置如下参数：
    1）设置测试次数
        TestNum=50;
    2）输入的图像以及输入模式
        img1='new2binarytreed-20.tif'; %输入图像I1
        img2='new2binarytreed-21.tif'; %输入图像I2
        DefTransType=2; %图像的输入模式。1表示用输入的M和V对I1做变换后得到I2，用I1和I2作为源图像进行运动估计；2表示I1和I2均直接读入。
    3）添加噪声
        DefNoisy=1; %在图像中添加噪声
        %噪声强度为Noise。由于添加噪声后图像的SNR与图像本身的内容相关，所以通过计算加躁后的SNR确定Noise的值。对于这组图像，SNR分别为0、5、10、15时的Noise取值已列在其后。下同。
        Noise=0.0018; % 0.07 0.02 0.005 0.0018
    4）设置算法
        algorithm=1;   %1表示使用我们的算法测试；2表示使用基于SIFT的算法；3表示使用[20]中的算法
        OR
        algorithm=2;

    运行Evaluation.m，得到相应的结果。输出的五组数据分别为MAE、MME、耗时、I1的SNR以及I2的SNR，每组数据的三个值分别为重复试验中得到的最小值、最大值和平均值。

    *********************************
    注：在不做批量测试而只对单幅图像做单次测试时，直接运行main.m或者SiftBased.m即可，上述2）和3）部分的参数可在main.m或者SiftBased.m中设置。下同。 
    *********************************
    
2.Forest 实验
    Evaluation.m中，使用%% nosie部分，注释%% 8 images部分，并设置如下参数：
    1）设置测试次数
        TestNum=50;
    2）输入的图像以及输入模式
        img1='forest3.bmp'; 
        DefTransType=1;
    3）添加噪声
        DefNoisy=1; 
        Noise=0.0016; % 0.07 0.019 0.005 0.0016
    4）设置仿射变换参数（单次试验时在main.m中设置）
        M=[0.05 0.01;0.01 0.06];
        V=[0.5 0.5];
    5）为减小在Deform过程中超出图像边界的像素点的损失，对图像做一些扩展（单次试验时在main.m中设置）
        DefExpand=1;
    6）设置算法
        algorithm=1;   
        OR
        algorithm=2;

    运行Evaluation.m，得到相应的结果。输出的五组数据分别为MAE、MME、耗时、I1的SNR以及I2的SNR，每组数据的三个值分别为重复试验中得到的最小值、最大值和平均值。

3.Lab 实验
    Evaluation.m中，使用%% nosie部分，注释%% 8 images部分，并设置如下参数：
    1）设置测试次数
        TestNum=50;
    2）输入的图像以及输入模式
        img1='lab.png'; 
        DefTransType=1;
    3）添加噪声
        DefNoisy=1; 
        Noise=0.0025; % 0.1 0.026 0.008 0.0025
    4）设置仿射变换参数
        M=[-0.01 -0.01;-0.03 0.02];
        V=[0.5 0.5];
    5）扩展图像
        DefExpand=1;
    6）设置算法
        algorithm=1; 

    运行Evaluation.m，得到相应的结果。输出的五组数据分别为MAE、MME、耗时、I1的SNR以及I2的SNR，每组数据的三个值分别为重复试验中得到的最小值、最大值和平均值。

4.对8幅图像在5组参数下的实验
    Evaluation.m中，使用%% 8 images部分，注释%% nosie部分，并设置如下参数：
    1）设置测试次数
        TestNum=8;
    2）输入模式
        DefTransType=1;
    4）设置仿射变换参数
        分别设置M和V为Table 1中的5组参数
    5）扩展图像
        DefExpand=1;
    6）设置算法
        algorithm=1; 
        OR
        algorithm=2;
        OR
        algorithm=3;

    运行Evaluation.m，得到相应的结果。输出的三组数据分别为MAE、MME、耗时，每组数据的8个值对应8幅图像。

    *********************************
    注：在这组实验中，由于图像以及变换参数变化范围较大，有时候需要调整一些main.m中的参数以获得更好的结果，如对于Par5，将计算层数由默认的5层提高为6层，对很多图像的结果有大幅提升。main.m中可能影响到计算结果的参数包括：
        MultiHeight %Multiscale计算层数
        DefBlurBefore %在对图像进行Deform前，是否进行模糊处理
        DefBlurAfter %在对图像进行Deform后，是否进行模糊处理，
        Blur %模糊处理的程度
    *********************************
