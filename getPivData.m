function pivData = getPivData(path)
    % 从速度场的txt文件中生成需要的pivData
%     path = './plane_v1.txt';
    data = importdata(path);
    pivData.x = data(:, 1)';
    pivData.y = data(:, 2)';
    pivData.ux = {};
    pivData.uy = {};
    for i = 1:30
        pivData.ux{i} = data(:, 3)';
        pivData.uy{i} = data(:, 4)';
    end
    pivData.type = '3D';
    pivData.domain = ~isnan(pivData.ux{1});