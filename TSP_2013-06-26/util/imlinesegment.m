function [data, btn] = imlinesegment()
hold on;
[x y btn] = ginput(1);

data = [];
i = 1;
while (btn == 1)
    data(end+1,1:2) = [x y];
    h(i) = plot(data(:,1), data(:,2), 'bx-', 'LineWidth', 3);
    drawnow;
    [x y btn] = ginput(1);
    i = i + 1;
end

for j=1:i-1
    delete(h(j))
end

end