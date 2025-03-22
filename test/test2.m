% 建立一個衛星模擬場景（假設有 satelliteScenario 函數）
scenario = satelliteScenario;

% 定義曼哈頓網格的行數和列數
numRows = 3;
numCols = 3;

% 衛星存放陣列（這裡假設有 satellite 函數可以依據給定位置創建衛星）
satellites = [];
for i = 1:numRows
    for j = 1:numCols
        % 設定衛星的軌道參數
        % 這裡以固定高度 500 km 為例，
        % 並透過調整緯度與經度來模擬曼哈頓式的排列（僅作示意，實際軌道設計會更複雜）
        lat = 40 + (i - 2) * 0.1;   % 緯度
        lon = -74 + (j - 2) * 0.1;  % 經度
        altitude = 500e3;          % 高度：500 km
        
        % 根據位置建立衛星（此函數名稱與參數僅為示例）
        sat = satellite(scenario, [lat, lon, altitude]);
        satellites = [satellites, sat];
    end
end

% 建立曼哈頓拓撲的連接
% 假設有 link 函數可以在兩個衛星之間建立連接
links = {};
for i = 1:numRows
    for j = 1:numCols
        idx = (i - 1) * numCols + j;
        % 若右邊有鄰近衛星，建立連線
        if j < numCols
            idxRight = idx + 1;
            lnk = link(satellites(idx), satellites(idxRight));
            links{end + 1} = lnk;
        end
        % 若下方有鄰近衛星，建立連線
        if i < numRows
            idxDown = idx + numCols;
            lnk = link(satellites(idx), satellites(idxDown));
            links{end + 1} = lnk;
        end
    end
end

% 執行模擬（根據工具箱可能需要額外設定）
simulate(scenario);

% 繪製或視覺化模擬場景
plot(scenario);
