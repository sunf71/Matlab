for i=1:1298
name = sprintf('d:\\sunf\\codes\\opencvtest\\iphone\\mb2s\\in%06d.jpg',i);
img = imread(name);
name = sprintf('d:\\sunf\\codes\\opencvtest\\iphone\\mb2sn\\in%06d.jpg',i);
imwrite(img,name);
end