function oim1_lab= rgb2lab(oim1)

cform = makecform('srgb2lab');
oim1_lab = lab2double(applycform(oim1,cform));
oim1_lab(:,:,1) = oim1_lab(:,:,1)*2.55;
oim1_lab(:,:,2) = oim1_lab(:,:,2) + 128;
oim1_lab(:,:,3) = oim1_lab(:,:,3) + 128;