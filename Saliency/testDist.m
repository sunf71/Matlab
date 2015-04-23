clear;

simg = imread('p1.png');

timg = imread('p2.png');

ssd1 = sum(sum(sum(abs(double(simg)-double(timg)))))/size(simg,1)/size(simg,2)/3/255;
mapping=getmapping(8,'u2');
H1=lbp(simg,1,8,mapping,'nh'); %LBP histogram in (8,1) neighborhood

%using uniform patterns
subplot(4,1,1),stem(H1);

H2=lbp(timg,1,8,mapping,'nh');
subplot(4,1,2),stem(H2);

lbp_dist1 = bhattacharyya(H1,H2);

simg = imread('p2.png');
timg = imread('p3.png');
ssd2 = sum(sum(sum(abs(double(simg)-double(timg)))))/size(simg,1)/size(simg,2)/3/255;

H1=lbp(simg,1,8,mapping,'nh'); %LBP histogram in (8,1) neighborhood
%using uniform patterns
subplot(4,1,3),stem(H1);

H2=lbp(timg,1,8,mapping,'nh');
subplot(4,1,4),stem(H2);

lbp_dist2 = bhattacharyya(H1,H2);

lbpd_ratio = lbp_dist2/lbp_dist1;
ssd_ratio = ssd2/ssd1;
fprintf('lbp dist ratio %f ssd dist ratio %f\n',lbpd_ratio,ssd_ratio);