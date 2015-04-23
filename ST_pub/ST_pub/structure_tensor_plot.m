function hST = structure_tensor_plot(Sdata)

% structure_tensor_plot - given a structure tensor, visualize as an ellipse %%%%
% 
% Example:
%   Sdata = partial_derivative_to_structure_tensor_form([1, 0.5, 0]);
%   structure_tensor_plot(Sdata);
% 
% Author: Shawn Arseneau
% Created: September 20, 2006
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    [e1,e2,e3,l1,l2,l3] = eigen_decomposition(Sdata);
    
    maxL = max([l1,l2,l3]);
    if maxL==0
        if nargout==1
            hST = [];
        end
       return; 
    else
        l1 = l1/maxL;
        l2 = l2/maxL;
        l3 = l3/maxL;
    end
    
    if l1>0 && l2==0 && l3==0  %--- single line
        [stheta, sphi, srho] = cart2sph(e3(1), e3(2), e3(3));
        [x, y, z] = sph2cart(stheta, sphi, 1); 
        plot3([-x, x], [-y, y], [-z, z], 'b', 'LineWidth', 3);
        hold on;
        maxDim = max([x,y,z]);
        axis equal;
        axis([-maxDim maxDim -maxDim maxDim -maxDim maxDim]);     
    else    
        [xe, ye, ze] = ellipsoid(0, 0, 0, l1, l2, l3, 50);          
        eData = [xe(:), ye(:), ze(:)];    
        eData = rotatePoints(e1, eData);
        hEll = surf(xe, ye, ze, 'FaceColor', 'b', 'EdgeColor', 'none', 'FaceAlpha', 0.85);
        hold on;
        axis equal;
        maxDim = max([1, eData(:)']);
        axis([-maxDim maxDim -maxDim maxDim -maxDim maxDim]);     
    end
    
    xlabel('X');
    ylabel('Y');
    zlabel('Z');
    grid on;
    if nargout==1
       hST = hEll; 
    end
    