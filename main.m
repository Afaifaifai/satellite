run('para.m')

sc = satelliteScenario(start_time, stop_time, sample_time, 'AutoSimulate', false);
restart(sc);

sat = walkerDelta(sc, semi_major_axis, inclination_deg, total_sat, P, f, Name="S");

gs1 = groundStation(sc, lat1, lon1, 'Name', bs_name1);
gs2 = groundStation(sc, lat2, lon2, 'Name', bs_name2);
uav1 = platform(sc, trajectory1, 'Name', uav_name1);
uav2 = platform(sc, trajectory2, 'Name', uav_name2);

wgs84 = referenceEllipsoid('wgs84');
[x, y, z] = geodetic2ecef(wgs84, lat1, lon1, bs_alt1);
ecefCoords = [x, y, z];
epoch = datetime('now');  % 確保使用 UTC 時間
utcArray = datevec(epoch);
gcrfCoords = ecef2eci(utcArray, ecefCoords);

sat1 = sat(1);
[pos, vel] = states(sat1);

% props = properties(sat(i));
% disp(props);
% disp(sat(i).OrbitPropagator);
% 設定模擬時間參數
endTime = sc.SimulationTime + hours(0.1);  % 模擬 1 小時後的時間
disp(sc.SimulationTime);
while sc.SimulationTime < endTime
    advance(sc);   % 推進模擬時間
    disp(sc.SimulationTime);
    disp(states(sat(1)))
    disp(1);
end

play(sc)