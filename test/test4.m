run('para.m')

sc = satelliteScenario(start_time, stop_time, sample_time, 'AutoSimulate', false);
restart(sc);

sat = walkerDelta(sc, semi_major_axis, inclination_deg, total_sat, P, f, Name="S");

gs1 = groundStation(sc, lat1, lon1, 'Name', bs_name1);
gs2 = groundStation(sc, lat2, lon2, 'Name', bs_name2);
% uav1 = platform(sc, trajectory1, 'Name', uav_name1);
% uav2 = platform(sc, trajectory2, 'Name', uav_name2);
uav1 = groundStation(sc, uav_lat1, uav_lon1, 'Altitude', uav_alt1, 'Name', uav_name1);
uav2 = groundStation(sc, uav_lat2, uav_lon2, 'Altitude', uav_alt2, 'Name', uav_name2);

wgs84 = referenceEllipsoid('wgs84');
[x, y, z] = geodetic2ecef(wgs84, lat1, lon1, bs_alt1);
ecefCoords = [x, y, z];
epoch = datetime('now');  % 確保使用 UTC 時間
utcArray = datevec(epoch);
gcrfCoords1 = ecef2eci(utcArray, ecefCoords);

[x, y, z] = geodetic2ecef(wgs84, lat2, lon2, bs_alt2);
ecefCoords = [x, y, z];
epoch = datetime('now');  % 確保使用 UTC 時間
utcArray = datevec(epoch);
gcrfCoords2 = ecef2eci(utcArray, ecefCoords);


r_mag = norm(gcrfCoords1);         % m
fprintf('地心距離：%.3f km\n', r_mag/1e3); % 證明地心到Coords1的距離等於地球半徑，因此GCRF的座標的1代表1公尺

display(gcrfCoords1);
display(gcrfCoords2);
display(gcrfCoords1 - gcrfCoords2);
% display(gcrfCoords1(1))
% display(gcrfCoords1(1) * 10)


%% 測試阻擋距離
% 1) 計算衛星間直線的最短距離到地心
p1 = gcrfCoords1;
p2 = gcrfCoords2;
v  = p2 - p1;                       % 方向向量
t0 = -dot(p1, v) / dot(v, v);      % 無限直線投影係數
% 若要針對「線段」，就夾限 t0 在 [0,1]：
t  = max(0, min(1, t0));           
closest = p1 + t * v;              % 線段上最接近地心的點
d_min = norm(closest);             % 地心到路徑的最短距離
fprintf('%.10f\n', d_min);
% 3) 判斷是否被阻擋
if d_min < block_radius
    disp('Warning: link blocked by Earth + atmosphere');
else
    disp('Link is clear');
end
%% 
sat1 = sat(1);
[pos, vel] = states(sat1);