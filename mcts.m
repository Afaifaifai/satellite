function bestPath = mcts(graph, startNode, endNode, iterations)
    % mctsShortestLoss 使用蒙地卡羅樹搜尋來找出從 startNode 到 endNode 的最小 Loss 路徑，
    % 且路徑長度不超過 maxHops+1（包含起始與終點）。
    %
    % 輸入：
    %   graph     - 鄰接矩陣，graph(i,j) 表示從 i 到 j 的 Loss，無連線者設定為 Inf
    %   startNode - 起點節點 (例如基地台)
    %   endNode   - 終點節點 (例如基地台)
    %   maxHops   - 允許的最大跳數（例如 3 表示路徑中有 4 個節點）
    %   iterations- MCTS 迭代次數
    %
    % 輸出：
    %   bestPath  - 估算出的最小 Loss 路徑 (節點編號序列)
    
    % 定義樹中節點的結構
    % 每個節點包含：state, path, cost, visits, totalCost, parent
    maxHops = Inf;
    node.state = startNode;        % 當前節點
    node.path = startNode;         % 從起點到目前為止的節點路徑
    node.cost = 0;                 % 累積 Loss
    node.visits = 0;               % 被訪問次數
    node.totalCost = 0;            % 累積模擬成本
    node.parent = 0;               % 父節點索引 (根節點的 parent 為 0)
    tree(1) = node;                % 初始化樹，第一個節點為根
    rootIdx = 1;
    
    bestCost = Inf;
    bestPath = [];
    
    for iter = 1:iterations
        %% 選擇 (Selection)
        currentIdx = rootIdx;
        while true
            currentNode = tree(currentIdx);
            % 若已到達終點（且非單一節點路徑）或已達到最大跳數則停止
            if (currentNode.state == endNode && length(currentNode.path) > 1) || (length(currentNode.path) == maxHops+1)
                break;
            end
            
            % 查看此節點是否有子節點（在樹中已展開的）
            childrenIdx = find([tree.parent] == currentIdx);
            if isempty(childrenIdx)
                break;
            else
                % 以 UCB（上置信界）方式選擇子節點，因為我們追求最小 Loss，
                % 這裡用平均成本減去探索項 (較低成本代表較好)
                bestValue = Inf;
                bestChild = [];
                for idx = childrenIdx
                    child = tree(idx);
                    if child.visits == 0
                        ucb = -Inf; % 尚未探索的節點先被選取
                    else
                        avgCost = child.totalCost / child.visits;
                        ucb = avgCost - 1.41 * sqrt(log(currentNode.visits+1) / child.visits);
                    end
                    if ucb < bestValue
                        bestValue = ucb;
                        bestChild = idx;
                    end
                end
                if isempty(bestChild)
                    break;
                else
                    currentIdx = bestChild;
                end
            end
        end
        
        %% 擴展 (Expansion)
        currentNode = tree(currentIdx);
        if ~(currentNode.state == endNode && length(currentNode.path) > 1) && (length(currentNode.path) < maxHops+1)
            % 找出從當前狀態可以移動到的下一個節點
            possibleMoves = find(~isinf(graph(currentNode.state, :)));
            % 可選擇性：避免循環（若不希望重複走過同一節點）
            possibleMoves = setdiff(possibleMoves, currentNode.path);
            if isempty(possibleMoves)
                newIdx = currentIdx; % 無法擴展，直接用當前節點
            else
                % 隨機選取一個動作擴展新節點
                nextState = possibleMoves(randi(length(possibleMoves)));
                newNode.state = nextState;
                newNode.path = [currentNode.path, nextState];
                newNode.cost = currentNode.cost + graph(currentNode.state, nextState);
                newNode.visits = 0;
                newNode.totalCost = 0;
                newNode.parent = currentIdx;
                tree(end+1) = newNode;
                newIdx = length(tree);
            end
        else
            newIdx = currentIdx;
        end
        
        %% 模擬 (Simulation)
        % 從新擴展的節點開始，隨機模擬直到到達終點或達到最大跳數
        simNode = tree(newIdx);
        while ~(simNode.state == endNode && length(simNode.path) > 1) && (length(simNode.path) < maxHops+1)
            possibleMoves = find(~isinf(graph(simNode.state, :)));
            possibleMoves = setdiff(possibleMoves, simNode.path);
            if isempty(possibleMoves)
                break;
            end
            nextState = possibleMoves(randi(length(possibleMoves)));
            simNode.path = [simNode.path, nextState];
            simNode.cost = simNode.cost + graph(simNode.state, nextState);
            simNode.state = nextState;
        end
        simCost = simNode.cost;
        
        % 若模擬路徑成功到達終點且成本較低，更新最佳解
        if simNode.state == endNode && simCost < bestCost
            bestCost = simCost;
            bestPath = simNode.path;
        end
        
        %% 回傳 (Backpropagation)
        % 從擴展節點向上更新沿途節點的統計資料
        idx = newIdx;
        while idx ~= 0
            tree(idx).visits = tree(idx).visits + 1;
            tree(idx).totalCost = tree(idx).totalCost + simCost;
            idx = tree(idx).parent;
        end
    end
    
    fprintf('Best path found: ');
    disp(bestPath);
    fprintf('Best total loss: %.2f\n', bestCost);
    
    end
    