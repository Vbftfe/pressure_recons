function coordinates = get_rect_cross_dots(vertexes)
% 获得由两个对角线顶点定义的矩形内的网格线与该对角线所有交点的坐标
% vertexes: 顶点坐标构成的数组 [v1_x, v1_y, v2_x, v2_y]
idx = 1;
if vertexes(1) > vertexes(3)
    temp = vertexes(1:2);
    vertexes(1:2) = vertexes(3:4);
    vertexes(3:4) = temp;
end

if vertexes(3) - vertexes(1) == 1
    coordinates = [];
    return;
end

for i = (vertexes(1)+1):(vertexes(3)-1)
    coord_y = (vertexes(4) - vertexes(2))/(vertexes(3) - vertexes(1))*(i - vertexes(1)) + vertexes(2);
    coordinates(idx, :) = [i, coord_y];
    idx = idx + 1;
end

for j = (vertexes(2)+1):(vertexes(4)-1)
    coord_x = (vertexes(3) - vertexes(1))/(vertexes(4) - vertexes(2))*(j - vertexes(2)) + vertexes(1);
    coordinates(idx, :) = [coord_x, j];
    idx = idx + 1;
end

coordinates = unique(coordinates, 'row');