%% 常量
mu = 3.986004418e14;                        % 地球引力參數 (m^3/s^2)
earth_radius = 6378137;                     % 地球赤道半徑 (m)
atm_height   = 500e3;                       % 大氣層高度 (m)
block_radius = earth_radius + atm_height;   % 地星到大氣層
h = 550e3;                                  % 衛星高度 (m)
semi_major_axis = earth_radius + h;         % 半長軸 (m) a
e = 0;                                      % 偏心率 (圓形軌道) e
inclination_deg = 53;                       % 傾角轉換為角度
inclination_rad = deg2rad(inclination_deg); % 傾角 (弧度)
omega = 0;                                  % 近地點幅角 (弧度)
omega_deg = rad2deg(omega);                 % 近地點幅角轉換為角度

start_time = datetime("yesterday", "TimeZone", "local");
stop_time = start_time + days(0.5);
sample_time = 60;

%% 星座參數
P = 24;                         % 軌道面數
S = 66;                         % 每面衛星數    
f = 11;                         % 相位因子
total_sat = P * S;              % 總衛星數

%% Base Station, UAV 參數
% 產生隨機經緯度（經度範圍 -180~180，緯度範圍 -90~90），假設高度固定為0公尺
% lat1 = -90 + 180*rand;
% lon1 = -180 + 360*rand;
% lat2 = -90 + 180*rand;
% lon2 = -180 + 360*rand;
lat1 = 23.45;
lon1 = -46.67;
lat2 = -23.45;
lon2 = 134.33;
bs_alt1 = 0;
bs_alt2 = 0;
bs_name1 = 'BS1';
bs_name2 = 'BS2';

% % 定義圓周上的點數
% nPoints = 10;
% angles = linspace(0, 2*pi, nPoints+1);
% angles(end) = [];  % 移除重複的起始點
% uav_radius = 10;  
% % 計算圓形軌跡上的航點
% uav_lats1 = lat1 + uav_radius * cos(angles);
% uav_lons1 = lon1 + uav_radius * sin(angles);
% uav_lats2 = lat2 + uav_radius * cos(angles);
% uav_lons2 = lon2 + uav_radius * sin(angles);
% % disp([uav_lats1(:), uav_lons1(:), repmat(uav_alt1, nPoints, 1)])
% trajectory1 = geoTrajectory([uav_lats1(:), uav_lons1(:), repmat(uav_alt1, nPoints, 1)], linspace(0, 3600, nPoints));
% trajectory2 = geoTrajectory([uav_lats2(:), uav_lons2(:), repmat(uav_alt2, nPoints, 1)], linspace(0, 3600, nPoints));

uav_lat1 = lat1;
uav_lon1 = lon1;
uav_lat2 = lat2;
uav_lon2 = lon2;
uav_alt1 = 1500;    % 高度 (m)
uav_alt2 = 1500;    % 高度 (m)
uav_name1 = 'UAV1 with BS1';
uav_name2 = 'UAV2 with BS2';


probs_cache_file_path = 'probs.mat';
sim_flag = false;
save_flag = true;