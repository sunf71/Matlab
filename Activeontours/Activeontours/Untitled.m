clear;
Start = 1;
End = 19;
sigma = 1.5;
epsilon = 1.5;
mu = 0.04;
lambda = 5;
alf = 1.5;
N = 300;
c0 = 4;
PlotRate = 40;
for i= Start : End;
    name = sprintf('warpErr%d.jpg',i);
    cimg = imread(name);
    simg = imresize(cimg,.25);
   
    LevelSetEvolutionWithoutReinitialization(simg,sigma,epsilon,mu,lambda,alf,c0,N,PlotRate,mask);
end