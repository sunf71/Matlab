function voxel_plot(pos, color, alpha);

%VOXEL function to draw a 3-D voxel of [1,1,1] in a 3-D plot %%%%%%%%%%%%%%
%
% Usage:
%   voxel(pos,color,alpha);
%
% Example:
%   voxel([2 3 4]);
%   axis([0 10 0 10 0 10]);
%

%   Created: Suresh Joel, Apr 15,2003, 
%   Updated  Suresh Joel, Feb 25, 2004
%   Updated: Shawn Arseneau, Sep.20, 2006
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if nargin<1 || nargin>3
       error('Incorrect number of inputs to voxel_plot'); 
    end
    if nargin<3
        alpha = 1;
    end
    if nargin<2
        color = [0.5, 0.5, 0.5];
    end
    d = [1,1,1];

x=[pos(1)+[0 0 0 0 d(1) d(1) d(1) d(1)]; ...
   pos(2)+[0 0 d(2) d(2) 0 0 d(2) d(2)]; ...
   pos(3)+[0 d(3) 0 d(3) 0 d(3) 0 d(3)]]';

for n=1:3
    if n==3
        x=sortrows(x,[n,1]);
    else
        x=sortrows(x,[n n+1]);
    end;
    temp=x(3,:);
    x(3,:)=x(4,:);
    x(4,:)=temp;
    h=patch(x(1:4,1),x(1:4,2),x(1:4,3),color);
    set(h,'FaceAlpha',alpha);
    set(h,'EdgeAlpha',0.05);
    temp=x(7,:);
    x(7,:)=x(8,:);
    x(8,:)=temp;
    h=patch(x(5:8,1),x(5:8,2),x(5:8,3),color);
    set(h,'FaceAlpha',alpha);
    set(h,'EdgeAlpha',0.05);
end;
