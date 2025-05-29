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

% start_time = datetime("yesterday", "TimeZone", "local");
start_time = datetime(2025, 5, 29, 17, 27, 0, 'TimeZone', 'local');

stop_time = start_time + days(1);
sample_time = 60 * 60 * 2; % seconds

%% 星座參數
P = 6;                         % 軌道面數
S = 11;                         % 每面衛星數    
f = 2;                         % 相位因子
total_sat = P * S;              % 總衛星數

lat1 = 23.45;
lon1 = -46.67;
lat2 = -23.45;
lon2 = 134.33;
bs_alt1 = 0;
bs_alt2 = 0;
bs_name1 = 'BS1';
bs_name2 = 'BS2';




sc = satelliteScenario(start_time, stop_time, sample_time, 'AutoSimulate', false);
restart(sc);

sats = walkerDelta(sc, semi_major_axis, inclination_deg, total_sat, P, f);
% sats = walkerDelta(sc, 7151000, 86.4, 66, 6, 2, Name="Iridium");
disp(length(sats));

gs1 = groundStation(sc, lat1, lon1, 'Name', bs_name1);
gs2 = groundStation(sc, lat2, lon2, 'Name', bs_name2);

uavs = walkerDelta(sc, earth_radius + 1e3, inclination_deg, 132, 12, 2, Name="UAV_Network");

epoch  = datetime(start_time);

% for i = 2:numel(sats)
%     disp(accessStatus(access(sats(1), sats(i)), start_time));
% end

% endTime = sc.SimulationTime + hours(0.1);  % 模擬 1 小時後的時間
endTime = sc.SimulationTime + days(1);
% end_Time = stop_time;
disp(sc.SimulationTime);


accesses = to_accesses(gs1, gs2, sats, uavs);
idx = 1;
while sc.SimulationTime < endTime
    coords = to_coords(gs1, gs2, sats, uavs, sc.SimulationTime);
    access_states = to_access_states(accesses, sc.SimulationTime);
    distances = coords2distances(coords);

    % 寫出 CSV 檔
    writematrix(coords,        sprintf('./data/iridium_coords_%03d.csv',        idx));
    writematrix(access_states, sprintf('./data/iridium_access_states_%03d.csv', idx));
    writematrix(distances,     sprintf('./data/iridium_distances_%03d.csv',     idx));

    disp(idx);
    idx = idx + 1;
    advance(sc);   % 推進模擬時間
end

play(sc)

% sat1 = sats(1);
% nowTime    = datetime('now');              % 取得目前系統時間

% for i = 2:numel(sats)
%     disp(i);
%     ac = access(sat1, sats(i));
%     disp(accessStatus(ac, start_time));
% end


function coords = to_coords(gs1, gs2, sats, uavs, time)
    % to_coords 將所有節點 (2 地面站 + sats + uavs) 在指定 time 的慣性座標取出
    function coord = transform_gcrf_coord(gs)
        wgs84 = referenceEllipsoid('wgs84');
        [x, y, z] = geodetic2ecef(wgs84, gs.Latitude, gs.Longitude, gs.Altitude);
        coord     = ecef2eci(datevec(time), [x, y, z]);  % 1×3
    end

    nSat = numel(sats);
    nUav = numel(uavs);
    N    = 2 + nSat + nUav;

    coords = zeros(N, 3);
    % 地面站
    coords(1, :) = transform_gcrf_coord(gs1);
    coords(2, :) = transform_gcrf_coord(gs2);

    % 衛星 sats
    for k = 1:nSat
        p = states(sats(k), time, "CoordinateFrame", "inertial");
        coords(2 + k, :) = p.';  % 第 3..(2+nSat)
    end

    % 無人機 uavs（假設也是 satellite 物件）
    for k = 1:nUav
        p = states(uavs(k), time, "CoordinateFrame", "inertial");
        coords(2 + nSat + k, :) = p.';  % 第 (3+nSat)..(2+nSat+nUav)
    end

    disp("Finishing transforming coordinates.");
end



function accesses = to_accesses(gs1, gs2, sats, uavs)
    % to_accesses 預先建立所有節點 pairwise 的 Access 物件
    nSat = numel(sats);
    nUav = numel(uavs);
    N    = 2 + nSat + nUav;  % 總節點數目

    accesses = cell(N, N);

    % 1↔2: 地面站互連
    accesses{1,2} = access(gs1, gs2);
    accesses{2,1} = accesses{1,2};

    % 地面站 ↔ sats
    for k = 1:nSat
        idx = 2 + k;
        accesses{1, idx} = access(gs1,    sats(k));
        accesses{idx, 1} = accesses{1, idx};
        accesses{2, idx} = access(gs2,    sats(k));
        accesses{idx, 2} = accesses{2, idx};
    end

    % 地面站 ↔ uavs
    for k = 1:nUav
        idx = 2 + nSat + k;
        accesses{1, idx} = access(gs1,    uavs(k));
        accesses{idx, 1} = accesses{1, idx};
        accesses{2, idx} = access(gs2,    uavs(k));
        accesses{idx, 2} = accesses{2, idx};
    end

    % sats ↔ sats
    for i = 1:nSat
        for j = i+1:nSat
            ii = 2 + i;
            jj = 2 + j;
            accesses{ii, jj} = access(sats(i), sats(j));
            accesses{jj, ii} = accesses{ii, jj};
        end
    end

    % uavs ↔ uavs
    for i = 1:nUav
        for j = i+1:nUav
            ii = 2 + nSat + i;
            jj = 2 + nSat + j;
            accesses{ii, jj} = access(uavs(i), uavs(j));
            accesses{jj, ii} = accesses{ii, jj};
        end
    end

    % sats ↔ uavs
    for i = 1:nSat
        for j = 1:nUav
            ii = 2 + i;
            jj = 2 + nSat + j;
            accesses{ii, jj} = access(sats(i), uavs(j));
            accesses{jj, ii} = accesses{ii, jj};
        end
    end

    disp("All Access objects created.");
end

function access_states = to_access_states(accesses, time)
    N = length(accesses);
    access_states = false(N, N);  % 初始化 N×N 的 logical 矩陣

    for i = 1:N
        for j = i+1:N
            s = accessStatus(accesses{i,j}, time);
            access_states(i,j) = s;
            access_states(j,i) = s;
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
