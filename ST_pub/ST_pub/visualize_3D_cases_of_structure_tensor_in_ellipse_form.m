function visualize_3D_cases_of_structure_tensor_in_ellipse_form(vizType)

% visualize_3D_cases_of_structure_tensor_in_ellipse_form - plots the surfel, curvel and ball form of the 3D structure tensor (demo) %%%%%%%%%%%%%%%%%%%%%
% 
% Example:
%  visualize_3D_cases_of_structure_tensor_in_ellipse_form('surfel');  %--- plots a surfel example
% 
%     output plots for a generic, 3D structure tensor, as well as the surfel, curvel and point versions
%     
%     Author: Shawn Arseneau
%     Created: May 16, 2006
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if nargin<1  
        vizType = 'gen_st';    
    end

    pos = [0,0,0];
    switch lower(vizType)
        case 'gen_st'  % --------------------------- Plot 3D Ellipse
            L = [15, 7, 2];
            aL = L(1)*1.5;

            [xe, ye, ze] = ellipsoid(0, 0, 0, L(1), L(2), L(3), 150);  
            surf(xe, ye, ze, 'FaceColor', 'interp', 'EdgeColor', 'none', 'FaceAlpha', 0.85);
            hold on;
            axis equal;
            camlight('headlight');
            lighting phong;
            colormap gray;
            arrow3D(pos, [0,0,aL], 'b');
            arrow3D(pos, [aL,0,0], 'r');
            arrow3D(pos, [0,aL,0], 'g');
            view(130,15);
            axis tight;    
            axis off;
            
        case 'surfel'  % --------------------------- Plot 3D Case for SURFEL
            hm = 40; incr = 10;
            [xgrid, ygrid] = meshgrid(-hm:incr:hm, -hm:incr:hm);

            tx = abs(xgrid)./max(xgrid(:));
            zgrid = zeros(size(xgrid))-tx.*0.5;
            surf(xgrid, ygrid, zgrid,'EdgeColor', [0.5 0.5 0.5], 'FaceColor', 'interp');
            hold on;
            [xe, ye, ze] = ellipsoid(0, 0, 0, 3, 6, 20, 150);  
            surf(xe, ye, ze, 'FaceColor', [0.25, 0.25 0.25], 'EdgeColor', 'none', 'FaceAlpha', 0.85);
            aL = 20*1.5;
            arrow3D(pos, [0,0,aL], 'b');
            arrow3D(pos, [aL,0,0], 'r');
            arrow3D(pos, [0,aL,0], 'g');
            axis equal; axis tight;
            view(135,12);
            camlight('headlight');    lighting phong;
            camlight right; camlight left;
            axis off;
            colormap autumn;
            
        case 'curvel' % ------------------------------- Plot 3D Case for CURVEL
            [xe, ye, ze] = ellipsoid(0, 0, 0, 25, 1, 1, 150);  
            surf(xe, ye, ze, 'EdgeColor', 'none');
            colormap hot;    
            [xe, ye, ze] = ellipsoid(0, 0, 0, 1, 8, 8, 150);  
            surf(xe, ye, ze, 'FaceColor', [0.85 0.85 0.85], 'EdgeColor', 'none', 'FaceAlpha', 0.7);
            hold on;    
            aL = 10;
            arrow3D(pos, [0,0,aL], 'b');
            arrow3D(pos, [aL,0,0], 'r');
            arrow3D(pos, [0,aL,0], 'g');
            axis([-12 12 -12 12 -10 10]);
            view(140,20);
            camlight headlight; lighting phong;
            axis off;
            colormap autumn;
            
        case 'ball'  % --------------------------------- Plot 3D Case for point (BALL)
            [x,y,z] = sphere(150); 
            radius = 7;
            x = x .* radius;
            y = y .* radius;
            z = z .* radius;

            surf(x,y,z,'FaceAlpha', 0.75, 'EdgeColor', 'none');
            hold on;    
            colormap gray;
            aL = 17;
            arrow3D(pos, [0,0,aL], 'b');
            arrow3D(pos, [aL,0,0], 'r');
            arrow3D(pos, [0,aL,0], 'g');
            axis equal;
            view(130,15);
            camlight left; lighting phong;
            axis off;
            
        otherwise
            tmsg = sprintf('otherwise,explanation_of_structure_tensor_in_ellipse_form:%s', lower(vizType)); 
            error(tmsg);
    end



    
    
    
    
    































