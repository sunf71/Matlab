function snr = SNR( I,In )
% 计算信号噪声比函数
% by Qulei
% I original signal
% In:noisy signal
% snr=10*log10(sigma2(I2)/sigma2(I2-I1))

[row,col,nchannel]=size(In);

snr=0;
if nchannel==1%gray image
% Ps=sum(sum((I-mean(mean(I))).^2));%signal power
% Pn=sum(sum((I-In).^2));%noise power
% snr=10*log10(Ps/Pn);

% Sn=abs(fft2(In)).^2;%noise power
% nA=sum(Sn(:))/numel(In);%noise average power
% Sf=abs(fft2(I)).^2;
% fA=sum(Sf(:))/numel(I);
% snr=10*log10(fA/nA);

In=reshape(I-In,1,row*col);
I=reshape(I,1,row*col);
Ps=std(I);%signal power
Pn=std(In);%noise power
snr=20*log10(Ps/Pn);

elseif nchannel==3%color image
for i=1:3
Ps=sum(sum((I(:,:,i)-mean(mean(I(:,:,i)))).^2));%signal power
Pn=sum(sum((I(:,:,i)-In(:,:,i)).^2));%noise power
snr=snr+10*log10(Ps/Pn);
end
snr=snr/3;
end

end

