function [x_pos, y_pos] = selectPoint(avgFrame,numberPoints,titleMessage)

h = figure(100); imshow(uint8(avgFrame)); title(titleMessage);
hold on
x_pos = zeros(1,numberPoints);
y_pos = zeros(1,numberPoints);
for k = 1:numberPoints
    [x_pos(k), y_pos(k)] = ginput(1);
    plot(x_pos(k),y_pos(k),'x','color','r','MarkerSize',5);
end
close(h)