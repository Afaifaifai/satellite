% function bestPath = mcts2(graph, coords, startNode, endNode, maxHops, iterations, probFunc)
%     % mcts_logprob 使用蒙地卡羅樹搜尋（Monte Carlo Tree Search, MCTS）
%     % 來最大化從 startNode 到 endNode 的邊機率乘積（product of probabilities），
%     % 且路徑中禁止重複經過節點（no cycles）。
%     %
%     % 輸入（Inputs）:
%     %   graph      - 鄰接矩陣（adjacency matrix），graph(i,j) 是節點 i 到 j 的距離（distance），
%     %                無連線者設為 Inf
%     %   startNode  - 起始節點編號（start node）
%     %   endNode    - 終點節點編號（end node）
%     %   maxHops    - 最大跳數（max hops），限制路徑長度 <= maxHops+1
%     %   iterations - MCTS 總迭代次數（number of iterations）
%     %   probFunc   - 函式句柄（function handle），將距離轉換為機率： p = probFunc(distance)
%     %
%     % 輸出（Output）:
%     %   bestPath   - 最佳路徑（node sequence）
    
%     % 節點（Node）結構定義
%     node.state    = startNode;        % 當前狀態（current state）
%     node.path     = startNode;        % 已走路徑（visited path）
%     node.visits   = 0;                % 拜訪次數（visit count）
%     node.totalLog = 0;                % 累積對數機率（sum of log-probabilities）
%     node.parent   = 0;                % 父節點索引（parent index）
%     tree(1)       = node;             % 樹初始化（initialize tree）
%     rootIdx       = 1;
    
%     bestLog = -Inf;    % 最佳平均 log-prob
%     bestPath = [];
    
%     % UCT 探索參數（exploration coefficient）
%     c = sqrt(2);
    
%     for iter = 1:iterations
%         %% 1. 選擇 (Selection)
%         currentIdx = rootIdx;
%         while true
%             curr = tree(currentIdx);
%             % 若已到終點或達最大跳數則停止選擇
%             if (curr.state == endNode && length(curr.path)>1) || length(curr.path)>=maxHops+1
%                 break;
%             end
%             % 找出已展開的子節點
%             children = find([tree.parent] == currentIdx); % 找出節點的parent是currentIdx的節點
%             if isempty(children)
%                 break;
%             end
%             % 用 UCT 公式選 child: UCT = avgLog + c*sqrt(log(N_parent)/N_i)
%             Np = curr.visits;
%             bestVal = -Inf;
%             bestChildIdx = children(1);
%             for idx = children
%                 child = tree(idx);
%                 if child.visits == 0
%                     uct = Inf;
%                 else
%                     avgLog = child.totalLog / child.visits;
%                     uct = avgLog + c * sqrt(log(Np)/child.visits);
%                 end
%                 if uct > bestVal
%                     bestVal = uct;
%                     bestChildIdx = idx;
%                 end
%             end
%             currentIdx = bestChildIdx;
%         end
    
%         %% 2. 擴展 (Expansion)
%         curr = tree(currentIdx);
%         if ~(curr.state==endNode && length(curr.path)>1) && length(curr.path)<maxHops+1
%             % 取得所有鄰居
%             % neigh = find(~isinf(graph(curr.state, :)));
%             % 過濾已訪問節點，禁止重複
%             % avail = setdiff(neigh, curr.path);

%             avail = s2s_(graph, coords, curr);
%             if ~isempty(avail)
%                 % 隨機選一個可行動作擴展
%                 nextState = avail(randi(length(avail)));
%                 newNode.state    = nextState; % 存的是原本矩陣的idx
%                 newNode.path     = [curr.path, nextState];
%                 % 計算對數機率增量 log p
%                 p = probFunc(graph(curr.state, nextState));
%                 newNode.visits   = 0;
%                 newNode.totalLog = 0;
%                 newNode.parent   = currentIdx;
%                 tree(end+1) = newNode;
%                 newIdx = length(tree);
%             else
%                 newIdx = currentIdx;
%             end
%         else
%             newIdx = currentIdx;
%         end
    
%         %% 3. 模擬 (Simulation)
%         sim = tree(newIdx);
%         simLog = 0;
%         while ~(sim.state==endNode && length(sim.path)>1) && length(sim.path)<maxHops+1
%             % neigh = find(~isinf(graph(sim.state, :)));
%             % avail = setdiff(neigh, sim.path);

%             avail = s2s_find_avail(graph, coords, sim);
%             if isempty(avail)
%                 break;
%             end
%             nextState = avail(randi(length(avail)));
%             % 累加對數機率
%             % p = probFunc(graph(sim.state, nextState));

%             p = calculate_probability();
%             simLog = simLog + log(p);
%             sim.path  = [sim.path, nextState];
%             sim.state = nextState;
%         end
    
%         %% 4. 回傳 (Backpropagation)
%         % 同時更新 visits 與 totalLog
%         idx = newIdx;
%         while idx~=0
%             tree(idx).visits   = tree(idx).visits + 1;
%             tree(idx).totalLog = tree(idx).totalLog + simLog;
%             idx = tree(idx).parent;
%         end
    
%         % 更新最佳路徑（最優平均 log-prob）
%         if sim.state == endNode
%             avgLog = tree(newIdx).totalLog / tree(newIdx).visits;
%             if avgLog > bestLog
%                 bestLog  = avgLog;
%                 bestPath = sim.path;
%             end
%         end
%     end
    
%     % 顯示結果
%     fprintf('Best path: ');
%     disp(bestPath);
% end


% function avail = s2s_find_avail(graph, coords, target_node) % sat to sat
%     neigh = find(~isinf(graph(target_node.state, :)));
%     diff_neigh = setdiff(neigh, target_node.path);

%     avail = [];
%     p1 = coords(target_node.state);
%     for idx = 1:numel(diff_neigh)
%         p2 = coords(diff_neigh(idx));
%         v  = p2 - p1;                       % 方向向量
%         t0 = -dot(p1, v) / dot(v, v);      % 無限直線投影係數
%         % 若要針對「線段」，就夾限 t0 在 [0,1]：
%         t  = max(0, min(1, t0));           
%         closest = p1 + t * v;              % 線段上最接近地心的點
%         d_min = norm(closest);             % 地心到路徑的最短距離
%         fprintf('%.10f\n', d_min);
%         % 3) 判斷是否被阻擋
%         if d_min < block_radius % link blocked by Earth + atmosphere
%             continue;
%         else
%             avail(end+1) = cur.state;
%         end
%     end
% end

% function p = calculate_probability() % poison
%     p = 1;
% end


function bestPath = mcts2(logP, startNode, endNode, iterations)
    % mcts_logprob 使用 MCTS 最大化起點到終點的機率乘積(以 log-p 加總形式)，
    % 禁止重複節點，直到到達終點或無路可走才停止。
    % 輸入:
    %   logP       - NxN 矩陣，logP(i,j) 為節點 i -> j 的 log(成功機率)，無連通為 Inf
    %   startNode  - 起始節點編號
    %   endNode    - 終點節點編號
    %   iterations - MCTS 模擬次數
    % 輸出:
    %   bestPath   - 最佳路徑節點編號序列
    
    % 直接跳到終點的小機率
    p_direct    = 0.01;
    
    % 初始化樹結構 (struct array)
    node.state    = startNode;
    node.path     = startNode;
    node.visits   = 0;
    node.totalLog = 0;
    node.parent   = 0;
    tree(1)       = node;
    rootIdx       = 1;
    
    bestLog  = -Inf;
    bestPath = [];
    
    % UCT 探索係數
    c = sqrt(2);
    
    for iter = 1:iterations
        fprintf('Iteration %d\n', iter);

        %% 1. Selection
        currentIdx = rootIdx;
        while true
            curr = tree(currentIdx);
            children = find([tree.parent]==currentIdx);
            if curr.state==endNode || isempty(children)
                break;
            end
            Np = curr.visits;
            bestVal = -Inf;
            for idx = children
                child = tree(idx);
                if child.visits == 0
                    uctVal = Inf;
                else
                    avgLog = child.totalLog / child.visits;
                    uctVal = avgLog + c * sqrt(log(Np)/child.visits);
                end
                if uctVal > bestVal
                    bestVal    = uctVal;
                    currentIdx = idx;
                end
            end
        end
    
        %% 2. Expansion
        curr = tree(currentIdx);
        if curr.state~=endNode
            neigh = find(~isinf(logP(curr.state,:)));
            avail = setdiff(neigh, curr.path);
            if ~isempty(avail)
                nextState = avail(randi(length(avail)));
                newNode.state    = nextState;
                newNode.path     = [curr.path, nextState];
                newNode.visits   = 0;
                newNode.totalLog = 0;
                newNode.parent   = currentIdx;
                tree(end+1)      = newNode;
                newIdx           = length(tree);
            else
                newIdx = currentIdx;
            end
        else
            newIdx = currentIdx;
        end
    
        %% 3. Simulation (roll-out)
        sim    = tree(newIdx);
        simLog = 0;
        while sim.state~=endNode
            % 嘗試直接跳到終點
            if rand() < p_direct && ~isinf(logP(sim.state,endNode))
                simLog = simLog + probP(endNode, sim.state);
                sim.state = endNode;
                sim.path(end+1) = endNode;
                break;
            end
            % 一般鄰居加權抽樣
            neigh = find(~isinf(logP(sim.state,:)));
            avail = setdiff(neigh, sim.path);
            if isempty(avail)
                simLog = -Inf;
                break;
            end
            % 計算線性機率權重
            pvec = exp(logP(sim.state, avail)); % 反指數化
            pvec = pvec / sum(pvec);
            cum = cumsum(pvec);
            sel = find(cum>=rand(),1);
            j   = avail(sel);
            % 累加 log-p
            simLog    = simLog + logP(sim.state, j);
            sim.state = j;
            sim.path(end+1) = j;
        end
    
        %% 4. Backpropagation
        idx = newIdx;
        while idx~=0
            tree(idx).visits   = tree(idx).visits + 1;
            tree(idx).totalLog = tree(idx).totalLog + simLog;
            idx = tree(idx).parent;
        end
    
        % 更新最佳路徑
        if sim.state==endNode
            avgLog = tree(newIdx).totalLog / tree(newIdx).visits;
            if avgLog > bestLog
                bestLog  = avgLog;
                bestPath = sim.path;
            end
        end
    end
    
    % 顯示結果
    fprintf('Best path: ');
    disp(bestPath);
    end
    
    
    
    