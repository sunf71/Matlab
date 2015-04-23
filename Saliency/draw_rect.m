function [state,result]=draw_rect(img,startPosition,windowSize,showOrNot)
% 函数调用：[state,result]=draw_rect(img,startPosition,windowSize,showOrNot)
% 函数功能：在图像画个长方形框
% 函数输入：img为原始的大图，可为灰度图，可为彩色图
%          startPosition 框的左上角在大图中的坐标(每行代表x,y坐标)，startPosition=[10,30],分别表示x,y为10,30
%          windowSize 框的大小 windowSize=[112,92] 分别表示宽、高
%          showOrNot 是否要显示结果?默认为显示出来?
% 函数输出：state -- 表示程序结果状态?
%          result - 结果图像数据 


if nargin < 4
    showOrNot = 1;
end

rgb = [255 0 0];                                 % 边框颜色
lineSize = 1;                                      % 边框大小，取1,2,3

windowSize(1,1)=windowSize(1,1);
windowSize(1,2) = windowSize(1,2);
if windowSize(1,2) > size(img,1) ||...
        windowSize(1,1) > size(img,2)
    state = -1;                                     % 说明窗口太大，图像太小，
    disp('the window size is larger then image...');
    return;
end

result = img;
if size(img,3) == 3
    for k=1:3
        for i=1:size(startPosition,1) %矩形框的总数
			if(startPosition(i,1)>=0 && startPosition(i,2)>=0)
				result(startPosition(i,2),startPosition(i,1):startPosition(i,1)+windowSize(i,1),k) = rgb(1,k); %画上边框  
				result(startPosition(i,2):startPosition(i,2)+windowSize(i,2),startPosition(i,1)+windowSize(i,1),k) = rgb(1,k);%画右边框
				result(startPosition(i,2)+windowSize(i,2),startPosition(i,1):startPosition(i,1)+windowSize(i,1),k) = rgb(1,k);  %画下边框   
				result(startPosition(i,2):startPosition(i,2)+windowSize(i,2),startPosition(i,1),k) = rgb(1,k);   %画左边框  
			
				if lineSize == 2 || lineSize == 3
					result(startPosition(i,2)+1,startPosition(i,1):startPosition(i,1)+windowSize(i,1),k) = rgb(1,k);  
					result(startPosition(i,2):startPosition(i,2)+windowSize(i,2),startPosition(i,1)+windowSize(i,1)-1,k) = rgb(1,k);
					result(startPosition(i,2)+windowSize(i,2)-1,startPosition(i,1):startPosition(i,1)+windowSize(i,1),k) = rgb(1,k);
					result(startPosition(i,2):startPosition(i,2)+windowSize(i,2),startPosition(i,1)-1,k) = rgb(1,k);
				
					if lineSize == 3
						result(startPosition(i,2)-1,startPosition(i,1):startPosition(i,1)+windowSize(i,1),k) = rgb(1,k);   
						result(startPosition(i,2):startPosition(i,2)+windowSize(i,2),startPosition(i,1)+windowSize(i,1)+1,k) = rgb(1,k);
						result(startPosition(i,2)+windowSize(i,2)+1,startPosition(i,1):startPosition(i,1)+windowSize(i,1),k) = rgb(1,k);
						result(startPosition(i,2):startPosition(i,2)+windowSize(i,2),startPosition(i,1)+1,k) = rgb(1,k);
					end
				end
			end
        end
    end
end

state = 1;

if showOrNot == 1
    figure;
	hold on;
    imshow(result);
end