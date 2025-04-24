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
            % disp(currentIdx);
            children = find([tree.parent]==currentIdx);
            % if (numel(children) >= 1)
            %     fprintf("  currIdx=%d, curr.state=%d, nc=%d %d\n", currentIdx, curr.state, numel(children), children(1));
            % end
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
                    uctVal = avgLog + c * sqrt(log(Np + 1)/child.visits);
                end
                if uctVal > bestVal
                    bestVal    = uctVal;
                    currentIdx = idx;
                end
            end
            % disp(currentIdx);
        end
        disp('select');
        %% 2. Expansion
        curr = tree(currentIdx);
        if curr.state~=endNode
            neigh = find(~isinf(logP(curr.state,:)));
            
            % disp(numel(neigh));
            avail = setdiff(neigh, curr.path);
            % disp(numel(avail));
            if ~isempty(avail)
                nextState = avail(randi(numel(avail)));
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
        % disp(tree);
    
        %% 3. Simulation (roll-out)
        sim    = tree(newIdx);
        simLog = 0;
        % disp(sim.state);

        while sim.state~=endNode
            % 嘗試直接跳到終點
            if rand() < p_direct && ~isinf(logP(sim.state,endNode))
                simLog = simLog + logP(endNode, sim.state);
                sim.state = endNode;
                sim.path(end+1) = endNode;
                break;
            end
            % 一般鄰居加權抽樣
            neigh = find(~isinf(logP(sim.state,:)));
            % neigh = find(~isinf(logP(curr.state,:)));
            % disp(numel(neigh));
            avail = setdiff(neigh, sim.path);
            % disp(numel(avail));
            % disp(numel(avail));
            if isempty(avail)
                simLog = -1e12;
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

    log_sum = 0;
    front = 1;
    for i = 2:numel(bestPath)
        log_sum = log_sum + logP(front, bestPath(i));
        front = bestPath(i);
    end
    disp(10^bestLog);
    disp(logP(2,130));
    end
    
    
    
    