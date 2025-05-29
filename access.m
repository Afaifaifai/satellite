%% 常量
mu = 3.986004418e14;                        % 地球引力參數 (m^3/s^2)
earth_radius = 6378137;                     % 地球赤道半徑 (m)
atm_height   = 500e3;                       % 大氣層高度 (m)
block_radius = earth_radius + atm_height;   % 地星到大氣層
h = 780e3;                                  % 衛星高度 (m)
semi_major_axis = earth_radius + h;         % 半長軸 (m) a
e = 0;                                      % 偏心率 (圓形軌道) e
inclination_deg = 86.4;                       % 傾角轉換為角度
inclination_rad = deg2rad(inclination_deg); % 傾角 (弧度)
omega = 0;                                  % 近地點幅角 (弧度)
omega_deg = rad2deg(omega);                 % 近地點幅角轉換為角度

start_time = datetime("yesterday", "TimeZone", "local");
stop_time = start_time + days(1);
sample_time = 60 * 60 * 2; % seconds

%% 星座參數
P = 6;                         % 軌道面數
S = 11;                         % 每面衛星數    
f = 2;                         % 相位因子
total_sat = P * S;              % 總衛星數

sc = satelliteScenario(start_time, stop_time, sample_time, 'AutoSimulate', false);
restart(sc);

sats = walkerDelta(sc, semi_major_axis, inclination_deg, total_sat, P, f);
% sats = walkerDelta(sc, 7151000, 86.4, 66, 6, 2, Name="Iridium");
disp(length(sats));

gs1 = groundStation(sc, lat1, lon1, 'Name', bs_name1);
gs2 = groundStation(sc, lat2, lon2, 'Name', bs_name2);

epoch  = datetime(start_time);

% for i = 2:numel(sats)
%     disp(accessStatus(access(sats(1), sats(i)), start_time));
% end

% endTime = sc.SimulationTime + hours(0.1);  % 模擬 1 小時後的時間
endTime = sc.SimulationTime + days(1);
% end_Time = stop_time;
disp(sc.SimulationTime);

idx = 1;
while sc.SimulationTime < endTime
    coords = to_coords(sats, gs1, gs2, sc.SimulationTime);
    access_states = to_access_states(sats, gs1, gs2, sc.SimulationTime);
    distances = coords2distances(coords);

    % 寫出 CSV 檔
    writematrix(coords,        sprintf('./data/iridium_coords_%03d.csv',        idx));
    writematrix(access_states, sprintf('./data/iridium_access_states_%03d.csv', idx));
    writematrix(distances,     sprintf('./data/iridium_distances_%03d.csv',     idx));

    disp(idx);
    idx = idx + 1;
    advance(sc);   % 推進模擬時間
end

% play(sc)

% sat1 = sats(1);
% nowTime    = datetime('now');              % 取得目前系統時間

% for i = 2:numel(sats)
%     disp(i);
%     ac = access(sat1, sats(i));
%     disp(accessStatus(ac, start_time));
% end


function coords = to_coords(sats, gs1, gs2, time)
    % epoch = datetime('now');  % 確保使用 UTC 時間
    function coord = transform_gcrf_coord(gs)
        wgs84 = referenceEllipsoid('wgs84');
        [x, y, z] = geodetic2ecef(wgs84, gs.Latitude, gs.Longitude, gs.Altitude);
        coord     = ecef2eci(datevec(time), [x, y, z]);  % 1×3
    end
    
    nSat   = numel(sats);
    coords = zeros(nSat + 2, 3);        % 每列一顆：2 ground + nSat

    coords(1,:) = transform_gcrf_coord(gs1);
    coords(2,:) = transform_gcrf_coord(gs2);

    for k = 1:nSat
        % states → 3×1×1 → squeeze → 3×1 → 轉列向量
        % time = datetime(2021,5,25,22,30,0);
        % display(k);
        p = states(sats(k), time, "CoordinateFrame","inertial");
        coords(k+2, :) = p.';            % 寫進第 k+2 列
    end

    disp("Finishing transforming code.");    
end


function access_states = to_access_states(sats, gs1, gs2, time)
    % sats: 1×n 衛星陣列
    % gs1, gs2: 地面站物件
    % time: datetime，要查詢的時間點

    n = numel(sats);
    N = n + 2;
    access_states = false(N, N);  % 初始化 N×N 的 logical 矩陣

    % 1↔2 (gs1 ↔ gs2)
    ac = access(gs1, gs2);
    access_states(1,2) = accessStatus(ac, time);
    access_states(2,1) = access_states(1,2);

    % 1↔(3:N) gs1 ↔ sats, 2↔(3:N) gs2 ↔ sats
    for i = 1:n
        idx = i + 2;
        % gs1 ↔ sats(i)
        ac = access(gs1, sats(i));
        access_states(1, idx) = accessStatus(ac, time);
        access_states(idx, 1) = access_states(1, idx);
        % gs2 ↔ sats(i)
        ac = access(gs2, sats(i));
        access_states(2, idx) = accessStatus(ac, time);
        access_states(idx, 2) = access_states(2, idx);
    end

    % sats ↔ sats
    for i = 1:n
        for j = 1:n
            if i ~= j
                ac = access(sats(i), sats(j));
                access_states(i+2, j+2) = accessStatus(ac, time);
            end
        end
    end

    % 強制對角線為 false
    access_states(1:N+1:end) = false;
end


function distances = coords2distances(coords)
%   coords2distances 計算節點間距離矩陣，利用對稱性減少計算量
%   coords : N×3 矩陣，每列為一個節點的 (x,y,z) 座標
%   distances      : N×N 距離矩陣，D(i,j)=D(j,i)，對角線為 0

    N = size(coords, 1);
    distances = zeros(N, N);         % 預先分配

    % 只計算 i<j 的部分
    for i = 1 : N-1
        for j = i+1 : N
            d = norm( coords(i,:) - coords(j,:) );
            distances(i,j) = d;
            distances(j,i) = d;      % 鏡射到另一半
        end
    end
end
