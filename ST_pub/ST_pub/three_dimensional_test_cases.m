function data3D = three_dimensional_test_cases(choice, maskSize)

% data3D - outputs one of several pre-defined 3D data sets %%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Example:
%   D1 = three_dimensional_test_cases('uniform');
% 
% Options:
%   choice = {'uniform','step_plane','sphere','line3d','random'}
%   maskSize = scalar denoting dimension. 3D size = maskSize x maskSize x maskSize
% 
% Author: Shawn Arseneau 
% Created: September 20, 2006
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if nargin<1  
       choice = 'uniform'; 
    end
    if nargin<2
        maskSize = 11;
    end
    midpt = ceil(maskSize/2);            

    switch(lower(choice))
        case 'uniform'
            data3D = zeros(maskSize, maskSize, maskSize);
            
        case 'step_plane'
            data3D = zeros(maskSize, maskSize, maskSize);
            data3D(1:midpt, :, :) = 1;
            
        case 'sphere'
            [xgrid, ygrid, zgrid] = meshgrid(-midpt:midpt, -midpt:midpt, -midpt:midpt);
            [stheta, sphi, rho] = cart2sph(xgrid, ygrid, zgrid);
            data3D = rho>=midpt;
            data3D = double(data3D);
            
        case 'line3d'
            data3D = ones(maskSize, maskSize, maskSize);
            data3D(midpt,midpt,:) = 0;
            
        case 'random'
            data3D = rand(maskSize, maskSize, maskSize);
            
        otherwise
            error('Unrecognized option in switch of three_dimensional_test_cases');
    end
    
%     %--- add some noise...
%     for d=1:maskSize        
%        data3D(:,:,d) = imnoise(data3D(:,:,d), 'gaussian');         
%     end
    
    



