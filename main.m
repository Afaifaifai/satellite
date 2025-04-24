% function main()
    run('para.m')

    
    % coords = to_coords(sats, gs1, gs2, datetime(start_time));

    if ~sim_flag && isfile(probs_cache_file_path)
        % S      = load(probs_cache_file_path, 'probs');
        % probs = S.probs;
        S      = load('coords.mat', 'coords');
        coords = S.coords;
        probs = to_probs(coords);
        fprintf("[INFO] Probabilities loaded from %s.\n", probs_cache_file_path);

    else
        sc = satelliteScenario(start_time, stop_time, sample_time, 'AutoSimulate', false);
        restart(sc);

        sats = walkerDelta(sc, semi_major_axis, inclination_deg, total_sat, P, f, Name="S");
        disp(length(sats));

        gs1 = groundStation(sc, lat1, lon1, 'Name', bs_name1);
        gs2 = groundStation(sc, lat2, lon2, 'Name', bs_name2);

        epoch  = datetime(start_time);
        coords = to_coords(sats, gs1, gs2, epoch);
        probs = to_probs(coords);
        % S      = load('coords.mat', 'coords');
        % coords = S.coords;
        % probs = to_probs(coords);
        fprintf("[INFO] Probabilities recomputed with to_probs().\n");

        if save_flag
            save(probs_cache_file_path, 'probs');
            fprintf("[INFO] Probabilities saved to %s.\n", probs_cache_file_path);
        end
    end

    % display(coords);

    probs(1,1) = Inf;
    probs(1,2) = Inf;
    probs(2,1) = Inf;
    probs(2,2) = Inf;
    % disp(probs);
    % for i = 3:1584
    %     neigh = find(~isinf(probs(i,:)));
    %     disp(numel(neigh));
    % end

    % disp(coords(10,:));

    % disp(probs(3, :));
    % disp(probs(3799,:));
    % disp(probs(2, :));

    mcts2(probs, 1, 2, 10);

    % disp(s2s_is_avail(coords(3), coords(7)));
    % disp(g2s_is_avail(coords(2), coords(3)));
    % for i = 1:numel(probs(1,:))
    %     if (probs(2, i) ~= Inf)
    %         % display(i);
    %     end
    % end


    % uav1 = platform(sc, trajectory1, 'Name', uav_name1);
    % uav2 = platform(sc, trajectory2, 'Name', uav_name2);

    % uav1 = groundStation(sc, uav_lat1, uav_lon1, 'Altitude', uav_alt1, 'Name', uav_name1);
    % uav2 = groundStation(sc, uav_lat2, uav_lon2, 'Altitude', uav_alt2, 'Name', uav_name2);
    % run('mcts2.m');

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


function probs = to_probs(coords)
    n_node = numel(coords(:,1));
    probs = zeros(n_node, n_node);

    for idx = 3:n_node
        if g2s_is_avail(coords(1,:), coords(idx,:))
            p = g2s_calculate_prob(coords(1), coords(idx));
            probs(1, idx) = log(p);
            probs(idx, 1) = log(p);
        else
            probs(1, idx) = Inf;
            probs(idx, 1) = Inf;
        end

        if g2s_is_avail(coords(2,:), coords(idx,:))
            p = g2s_calculate_prob(coords(2,:), coords(idx,:));
            probs(2, idx) = log(p);
            probs(idx, 2) = log(p);
        else
            probs(2, idx) = Inf;
            probs(idx, 2) = Inf;
        end

        for jdx = 3:n_node
            if (idx == jdx)
                probs(idx, jdx) = Inf;
                continue;
            end
            if s2s_is_avail(coords(idx,:), coords(jdx,:))
                p = s2s_calculate_prob(coords(idx,:), coords(jdx,:));
                probs(jdx, idx) = log(p);
                probs(idx, jdx) = log(p);
            else
                probs(jdx, idx) = Inf;
                probs(idx, jdx) = Inf;
            end
        end
    end
end

function avail = g2s_is_avail(gs_coord, sat_coord)
    earth_radius = 6378137;                     % 地球赤道半徑 (m)
    h = 550e3;

    r_gs  = norm(gs_coord);
    r_sat = norm(sat_coord);
    cos_theta = dot(gs_coord, sat_coord) / (r_gs * r_sat);
    % fprintf('%.4f %.4f', r_gs, r_sat);
    if cos_theta > (earth_radius / (earth_radius + h))
        avail = true;
    else
        avail = false;
    end
end

function avail = s2s_is_avail(sat1_coord, sat2_coord)
    earth_radius = 6378137;                     % 地球赤道半徑 (m)
    atm_height   = 200e3;                       % 大氣層高度 (m)
    block_radius = earth_radius + atm_height;

    v  = sat2_coord - sat1_coord;                       % 方向向量
    t0 = -dot(sat1_coord, v) / dot(v, v);      % 無限直線投影係數
    % 若要針對「線段」，就夾限 t0 在 [0,1]：
    t  = max(0, min(1, t0));           
    closest = sat1_coord + t * v;              % 線段上最接近地心的點
    d_min = norm(closest);             % 地心到路徑的最短距離
    % 3) 判斷是否被阻擋
    if d_min < block_radius % link blocked by Earth + atmosphere
        avail = false;
    else
        avail = true;
    end
end

function p = g2s_calculate_prob(gs_coord, sat_coord)
    p = rand;
end

function p = s2s_calculate_prob(sat_coord1, sat_coord2)
    p = rand;
end

% props = properties(sat(i));
% disp(props);
% disp(sat(i).OrbitPropagator);
% 設定模擬時間參數

% endTime = sc.SimulationTime + hours(0.1);  % 模擬 1 小時後的時間
% disp(sc.SimulationTime);
% while sc.SimulationTime < endTime
%     advance(sc);   % 推進模擬時間
%     disp(sc.SimulationTime);
%     disp(states(sat(1)))
%     disp(1);
% end

% play(sc)