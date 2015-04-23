function test_image(choice)

% test_image - outputs a rampStep, stepEdge, circle and uniform image for the demoStructureTensor.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     test_image(choice)
% 
%     Given 'choice', output simulated image
%     
% Example:
%  test_image(1);  %----- creates and saves a ramp-step-edge image
% 
%     Author: Shawn Arseneau
%     Created: Sept.15, 2006
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if nargin~=1   
        choice = 1;   
    end
    
    switch(choice)
        case 1 
            msize = 21;
            hmask = floor(msize/2);
            img = zeros(msize, msize);
            B2W = (0:1/(msize-1):1);
            for i=1:hmask
               img(i,:) = B2W; 
            end            
            imwrite(img, 'rampStep.pgm');
        
        case 2
            msize = 21;
            hmask = floor(msize/2);
            img = zeros(msize, msize);
            
            for i=1:hmask
               img(i,:) = 1; 
            end            
            imwrite(img, 'stepEdge.pgm');
            
            
        case 3   %---- black circle
            msize = 151;
            hmask = floor(msize/2);
            [xG, yG] = meshgrid(-hmask:hmask, -hmask:hmask);
            [thetaG, rhoG] = cart2pol(xG, yG);
            img = rhoG>hmask;
            imwrite(double(img), 'circle.pgm');
            
        case 4  %---- monochromatic
            msize = 61;
            img = ones(msize, msize).*(0.5);
            imwrite(double(img), 'oneColor.pgm');            
            
        otherwise
            error('otherwise reached in switch for test_image');        
    end







