function visualize_3D_data(data)

% visualize_3D_data - given an NxMxP data set, output as a collection of voxels
% 
% Example:
%    R = rand(5,4,3);
%    visualize_3D_data(R);
% 
% 
% Author: Shawn Arseneau
% Created: September 20, 2006
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if nargin~=1
       error('no data sent to visualize_3D_data'); 
    end
    if ndims(data)~=3
       error('data sent to visualize_3D_data not in NxMxP form'); 
    end
    
    [N,M,P] = size(data);    
    hold on;
    maxColor = max([1, max(data(:))]);
    colorData = ((data./maxColor).*0.5) + 0.5 ; %---- 0.5 ... 1
        
    for r=1:N
       for c=1:M 
           for d=1:P               
               vColor = [1,1,1].*colorData(c,r,d);
               voxel_plot([c,r,d],vColor, 1-data(c,r,d));
           end
       end
    end
    axis equal;
    axis tight;
    view(30,30);
    xlabel('X');
    ylabel('Y');
    zlabel('Z');
    
    
    
    
    
    




