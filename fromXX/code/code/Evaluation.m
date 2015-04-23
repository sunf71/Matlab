clc;
clear all;

TestNum=25;
algorithm=1; %1:ours; 2:Sift-based 3:PMI

   img1='new2binarytreed-20.tif';
   img2='new2binarytreed-21.tif';
%   img1='forest3.bmp';
%   img1='office.bmp';
%   
  DefTransType=2;
  
  DefExpand=0;
  DefNoisy=1; Noise=0.0018;
% 
  M=[0.05   0.010
     0.010   0.06];
  V=[0.5   0.5];
%   M=[-0.01   -0.01
%     -0.03   0.02];
%   V=[0.5   0.5];

%   M=[0.02   0.010
%      0.20   0.04];
%   V=[0.1   0.2];
  
%   M=[0.006 0.068
%     -0.034 0.012];
% V=[0.1 0.2];
% 
%  M=[0.1   0.120
%      0.090   -0.2];
%   V=[3   0.5];
% 
%   M=[0.07   -0.12
%      0.09   0.1];
%   V=[1   2];

%  M=[-0.3   0.120
%     0.120   -0.5];
%   V=[10.0   20.0];
  
EVA=1;
angs=zeros(1,TestNum);
mags=zeros(1,TestNum);
times=zeros(1,TestNum);
SNR1s=zeros(1,TestNum);
SNR2s=zeros(1,TestNum);

%% nosie
 
 TransType=DefTransType;
for testtime=1:TestNum
    DefTransType=TransType;
    
    if algorithm==1
        main;
    elseif algorithm==2
        SiftBased;
    end
    
    angs(testtime)=ang;
    mags(testtime)=mag;
    times(testtime)=time;
    if exist('snr1','var') && exist('snr2','var')
    SNR1s(testtime)=snr1;
    SNR2s(testtime)=snr2;
    end
    
    k=round(testtime/TestNum*100);
    if mod(k,2)==0
        clc;
        disp(['Complete ',num2str(k),'%']);
    end
end

clear EVA
clc;
disp([num2str(TestNum),' times of evalutaion completed']);
disp(['angs: ',num2str(min(angs)),' ',num2str(max(angs)),' ',num2str(mean(angs))]);
disp(['mags: ',num2str(min(mags)),' ',num2str(max(mags)),' ',num2str(mean(mags))]);
disp(['times: ',num2str(min(times)),' ',num2str(max(times)),' ',num2str(mean(times))]);
disp(['SNR1s: ',num2str(min(SNR1s)),' ',num2str(max(SNR1s)),' ',num2str(mean(SNR1s))]);
disp(['SNR2s: ',num2str(min(SNR2s)),' ',num2str(max(SNR2s)),' ',num2str(mean(SNR2s))]);

%% 8 images

%  TransType=DefTransType;
% for testtime=1:8
%     DefTransType=TransType;
%     img1=['8images\' num2str(testtime) '.jpg'];
%     if algorithm==1
%        main;
%     elseif algorithm==2
%        SiftBased;
%     elseif algorithm==3
%         radonPMI;
%     end
%     angs(testtime)=ang;
%     mags(testtime)=mag;
%     times(testtime)=time;
%     
%     clc;
%     disp(['Complete ',num2str(testtime),' pics']);
% end
% clear EVA
% clc;
% angs
% mags
% times