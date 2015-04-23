function [e1,e2,e3,l1,l2,l3] = eigen_decomposition(M)

% eigen_decomposition = uses eigs.m and appropriately labels output (for both 2D and 3D)  %%%%%%%%%%%%%%
%
%   [e1,e2,l1,l2] = eigen_decomposition(M)
%
%   Perform eigen-based decomposition of the input tensor M.
%   
%   INPUT:
%    M = structure tensor of the form [dxdx dxdy; dxdy dydy]
%
% Example:
% [e1,e2,l1,l2] = eigen_decomposition([1,0; 0,1]);
% 
%   Author: Shawn Arseneau
%   Created: August 22, 2005
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if (nargout~=4 && nargout~=6)
       error('input to eigen_decomposition must be a 2x2 matrix and there must be four output variables'); 
    end

    if ndims(M)==2
        [V D] = eigs(M);
        D = abs(D);
        l1 = D(1,1);
        l2 = D(2,2);
        e1 = V(:,1);
        e2 = V(:,2);

        if nargout==6
           l3 = D(3,3); 
           e3 = V(:,3);
           return;
        else %------- form [e1,e2,l1,l2]
            e3 = l1;
            l1 = l2;
        end
    elseif ndims(M)==4
        [rows, cols, st1, st2] = size(M);
        tMat = zeros(st1, st2);
        
        if nargout==6 %---- 3D info
            for r=1:rows
                for c=1:cols
                    tMat(:,:) = M(r,c,:,:);
                    [V D] = eigs(tMat);
                    D = abs(D);
                    l1(r,c) = D(1,1);
                    l2(r,c) = D(2,2);
                    l3(r,c) = D(3,3); 
                    e1(r,c,:) = V(:,1);
                    e2(r,c,:) = V(:,2);
                    e3(r,c,:) = V(:,3);
                end
            end
        else  %---- 2D info
            for r=1:rows
                for c=1:cols
                    tMat(:,:) = M(r,c,:,:);
                    [V D] = eigs(tMat);
                    D = abs(D);
                    e1(r,c,:) = V(:,1);
                    e2(r,c,:) = V(:,2);
                    e3(r,c) = D(1,1);
                    l1(r,c) = D(2,2);
                end
            end
        end
    else
        error('Input to eigen decomposition in wrong array amount (NxN) or (NxMx2x2)');
    end
    
    
    
    
    
    
    
    
    
    
    