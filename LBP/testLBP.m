       I=imread('in000244.jpg');
       mapping=getmapping(8,'u2'); 
       H1=LBP(I,1,8,mapping,'h'); %LBP histogram in (8,1) neighborhood
                                  %using uniform patterns
       subplot(2,1,1),stem(H1);

       H2=LBP(I);
       subplot(2,1,2),stem(H2);

       SP=[-1 -1; -1 0; -1 1; 0 -1; -0 1; 1 -1; 1 0; 1 1];
       I2=LBP(I,SP,0,'i'); %LBP code image using sampling points in SP
                           %and no mapping. Now H2 is equal to histogram
                           %of I2.