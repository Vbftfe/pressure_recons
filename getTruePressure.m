function truePressure = getTruePressure(path)
% path = './pressure.txt';
    data = importdata(path);
    truePressure = data(:, 3)';

