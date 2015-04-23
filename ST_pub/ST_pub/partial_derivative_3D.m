function [Ix,Iy,Iz] = partial_derivative_3D(data)

% partial_derivative_3D - average of absolute, successive difference along axes %%%
% 
% Given data of the form NxMxP, output the partial derivative information.
% This is calculated as the average neighbor absolute difference.
% 
% Example:
%   R = rand(5,3,4);
%   [Rx, Ry, Rz] = partial_derivative_3D(R);
% 
% Author: Shawn Arseneau
% Created: September 20, 2006
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if ndims(data)~=3
       error('data is not in NxMxP form'); 
    end
    if nargout~=3
       error('incorrect number of output parameters'); 
    end
    
    %--- difference along rows
    A = data(1:end-1, :, :);
    B = data(2:end, :, :);
    Iy = abs(A-B);
    numPoints = numel(Iy);
    Iy = sum(Iy(:))/numPoints;
    
    %--- difference along columns
    A = data(:, 1:end-1, :);
    B = data(:, 2:end, :);
    Ix = abs(A-B);
    numPoints = numel(Ix);
    Ix = sum(Ix(:))/numPoints;
    
    %--- difference along depth
    A = data(:, :, 1:end-1);
    B = data(:, :, 2:end);
    Iz = abs(A-B);
    numPoints = numel(Iz);
    Iz = sum(Iz(:))/numPoints;
    
    
    
















