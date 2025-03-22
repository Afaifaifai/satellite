% #grok
para;

% 計算平均運動量
n = sqrt(mu / a^3); % 平均運動量 (rad/s)

% 計算每個衛星的RAAN和初始平均近點角
Omega_all = zeros(total_sat, 1);
M0_all = zeros(total_sat, 1);

idx = 1;
for k = 0:P-1
    Omega_k = k * 2 * pi / P; % RAAN for plane k (radians)
    for j = 0:satellitesPerPlane-1
        % 計算平均近點角
        M0_kj = j * (2 * pi / satellitesPerPlane) + k * (f * 2 * pi / total_sat);
        M0_kj = mod(M0_kj, 2 * pi); % 確保在[0, 2*pi)範圍內
        
        Omega_all(idx) = Omega_k;
        M0_all(idx) = M0_kj;
        
        idx = idx + 1;
    end
end

% 計算單個衛星在時間t的ECI位置的函數
function r_ECI = get_position(m, t, a, e, i_rad, Omega_all, M0_all, n)
    M_t = mod(M0_all(m) + n * t, 2 * pi);
    v = M_t; % 對於e=0，真近點角等於平均近點角
    r_orb = a * [cos(v); sin(v); 0];
    
    cos_Omega = cos(Omega_all(m));
    sin_Omega = sin(Omega_all(m));
    cos_i = cos(i_rad);
    sin_i = sin(i_rad);
    
    R3_Omega = [cos_Omega, -sin_Omega, 0; sin_Omega, cos_Omega, 0; 0, 0, 1];
    R1_i = [1, 0, 0; 0, cos_i, -sin_i; 0, sin_i, cos_i];
    
    R = R3_Omega * R1_i;
    r_ECI = R * r_orb;
end

% 範例：計算t=0時所有衛星的位置
positions = zeros(3, total_sat);
for m = 1:total_sat
    positions(:, m) = get_position(m, 0, a, e, i_rad, Omega_all, M0_all, n);
end

disp('前5顆衛星在t=0的ECI位置 (單位：米)：');
disp(positions(:, 1:5));

% 模擬一天的運動，步長1分鐘
t_start = 0;
t_end = 86400; % 1天 (秒)
dt = 60; % 1分鐘
times = t_start:dt:t_end;

num_times = length(times);
positions_over_time = zeros(3, total_sat, num_times);

for ti = 1:num_times
    t = times(ti);
    for m = 1:total_sat
        positions_over_time(:, m, ti) = get_position(m, t, a, e, i_rad, Omega_all, M0_all, n);
    end
end