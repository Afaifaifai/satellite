function bestPath = mcts_logprob(graph, startNode, endNode, maxHops, iterations, probFunc)
    % mcts_logprob 使用蒙地卡羅樹搜尋（Monte Carlo Tree Search, MCTS）
    % 來最大化從 startNode 到 endNode 的邊機率乘積（product of probabilities），
    % 且路徑中禁止重複經過節點（no cycles）。
    %
    % 輸入（Inputs）:
    %   graph      - 鄰接矩陣（adjacency matrix），graph(i,j) 是節點 i 到 j 的距離（distance），
    %                無連線者設為 Inf
    %   startNode  - 起始節點編號（start node）
    %   endNode    - 終點節點編號（end node）
    %   maxHops    - 最大跳數（max hops），限制路徑長度 <= maxHops+1
    %   iterations - MCTS 總迭代次數（number of iterations）
    %   probFunc   - 函式句柄（function handle），將距離轉換為機率： p = probFunc(distance)
    %
    % 輸出（Output）:
    %   bestPath   - 最佳路徑（node sequence）
    
    % 節點（Node）結構定義
    node.state    = startNode;        % 當前狀態（current state）
    node.path     = startNode;        % 已走路徑（visited path）
    node.visits   = 0;                % 拜訪次數（visit count）
    node.totalLog = 0;                % 累積對數機率（sum of log-probabilities）
    node.parent   = 0;                % 父節點索引（parent index）
    tree(1)       = node;             % 樹初始化（initialize tree）
    rootIdx       = 1;
    
    bestLog = -Inf;    % 最佳平均 log-prob
    bestPath = [];
    
    % UCT 探索參數（exploration coefficient）
    c = sqrt(2);
    
    for iter = 1:iterations
        %% 1. 選擇 (Selection)
        currentIdx = rootIdx;
        while true
            curr = tree(currentIdx);
            % 若已到終點或達最大跳數則停止選擇
            if (curr.state == endNode && length(curr.path)>1) || length(curr.path)>=maxHops+1
                break;
            end
            % 找出已展開的子節點
            children = find([tree.parent] == currentIdx); % 找出節點的parent是currentIdx的節點
            if isempty(children)
                break;
            end
            % 用 UCT 公式選 child: UCT = avgLog + c*sqrt(log(N_parent)/N_i)
            Np = curr.visits;
            bestVal = -Inf;
            bestChildIdx = children(1);
            for idx = children
                child = tree(idx);
                if child.visits == 0
                    uct = Inf;
                else
                    avgLog = child.totalLog / child.visits;
                    uct = avgLog + c * sqrt(log(Np)/child.visits);
                end
                if uct > bestVal
                    bestVal = uct;
                    bestChildIdx = idx;
                end
            end
            currentIdx = bestChildIdx;
        end
    
        %% 2. 擴展 (Expansion)
        curr = tree(currentIdx);
        if ~(curr.state==endNode && length(curr.path)>1) && length(curr.path)<maxHops+1
            % 取得所有鄰居
            neigh = find(~isinf(graph(curr.state, :)));
            % 過濾已訪問節點，禁止重複
            avail = setdiff(neigh, curr.path);
            if ~isempty(avail)
                % 隨機選一個可行動作擴展
                nextState = avail(randi(length(avail)));
                newNode.state    = nextState;
                newNode.path     = [curr.path, nextState];
                % 計算對數機率增量 log p
                p = probFunc(graph(curr.state, nextState));
                newNode.visits   = 0;
                newNode.totalLog = 0;
                newNode.parent   = currentIdx;
                tree(end+1) = newNode;
                newIdx = length(tree);
            else
                newIdx = currentIdx;
            end
        else
            newIdx = currentIdx;
        end
    
        %% 3. 模擬 (Simulation)
        sim = tree(newIdx);
        simLog = 0;
        while ~(sim.state==endNode && length(sim.path)>1) && length(sim.path)<maxHops+1
            neigh = find(~isinf(graph(sim.state, :)));
            avail = setdiff(neigh, sim.path);
            if isempty(avail)
                break;
            end
            nextState = avail(randi(length(avail)));
            % 累加對數機率
            p = probFunc(graph(sim.state, nextState));
            simLog = simLog + log(p);
            sim.path  = [sim.path, nextState];
            sim.state = nextState;
        end
    
        %% 4. 回傳 (Backpropagation)
        % 同時更新 visits 與 totalLog
        idx = newIdx;
        while idx~=0
            tree(idx).visits   = tree(idx).visits + 1;
            tree(idx).totalLog = tree(idx).totalLog + simLog;
            idx = tree(idx).parent;
        end
    
        % 更新最佳路徑（最優平均 log-prob）
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
    