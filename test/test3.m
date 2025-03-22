run('para.m');  % 載入基本參數

%% 建立衛星模擬場景
scenario = satelliteScenario;

%% 建立 Walker Delta 衛星系統，使用 24x66 的 cell array 儲存衛星物件
satellites = cell(P, S);  % 建立一個 24x66 的 cell array

for plane = 0:(P-1)
    % 計算每個軌道面的升交點赤經 (RAAN)，均勻分佈於 360°
    RAAN = plane * (360 / P);
    
    for satInPlane = 0:(S-1)
        % 根據 Walker Delta 模式計算真近點角 (True Anomaly)
        % 公式： trueAnomaly = (360/S) * (satInPlane + plane*(f/P))
        trueAnomaly = mod((360 / S) * (satInPlane + plane*(f/P)), 360);
        
        % 定義衛星名稱
        satName = sprintf('Sat_%02d_%02d', plane+1, satInPlane+1);
        
        tleFile = "eccentricOrbitSatellite.tle";
        % 建立衛星物件（僅傳入名稱，不在此設定軌道）
        satObj = satellite(scenario, tleFile, 'Name', satName);
        
        % 組成 Keplerian 軌道參數向量 (確保是 1x6 的數值列向量)
        orbitParams = [a, e, i_deg, RAAN, omega_deg, trueAnomaly];
        orbitParams = orbitParams(:).';  % 確保是 1×6
        
        % 現在再設定軌道參數
        satObj.Orbit.Keplerian = orbitParams;
        
        % 將此衛星物件儲存到 cell array 的對應位置中
        satellites{plane+1, satInPlane+1} = satObj;
    end
end

%% 視覺化
viewer = satelliteScenarioViewer(scenario);

% 為每顆衛星新增地面軌跡（Ground Track）
for plane = 1:P
    for satInPlane = 1:S
        GroundTrack(satellites{plane, satInPlane});
    end
end
