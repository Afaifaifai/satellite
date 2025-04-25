function bestPath = mcts2(lambdaMat, startNode, endNode, iterations)
    %-----------------------------------------------------------------
    % MCTS (Poisson‑edge) with percentile ribbon & data logging
    % lambdaMat(i,j) = Poisson λ of edge i→j ;  Inf  = no link
    %-----------------------------------------------------------------
    
    p_direct  = 0.01;              % chance to jump straight to end
    c         = sqrt(2);           % UCT exploration coefficient
    NEG_INF   = -1e12;             % penalty for dead end (finite)
    
    %%--- root node struct ------------------------------------------
    node.state  = startNode;
    node.path   = startNode;
    node.depth  = 0;               % hop index (root =0)
    node.visits = 0;
    node.totalLg= 0;
    node.parent = 0;
    node.isFull = false;
    
    tree(1) = node;  rootIdx = 1;
    
    bestLg   = NEG_INF;  bestPath = [];
    
    %%--- logging vectors -------------------------------------------
    log_record  = NEG_INF*ones(1,iterations);   % roll‑out logP
    prob_record = zeros(1,iterations);          % roll‑out product prob
    step_record = zeros(1,iterations);          % roll‑out steps
    best_record = NEG_INF*ones(1,iterations);   % best log so far
    size_record = zeros(1,iterations);          % tree size
    
    %%================== MAIN MCTS LOOP =============================
    for it = 1:iterations
        %% 1‑Selection
        curIdx = rootIdx;
        while true
            cur = tree(curIdx);
            kids = find([tree.parent]==curIdx & ~[tree.isFull]);
            if cur.state==endNode || isempty(kids), break; end
            Np = cur.visits + 1; bestVal = -Inf; bestKid = kids(1);
            for k = kids
                ch = tree(k);
                if ch.visits==0, u = Inf;
                else, avg = ch.totalLg/ch.visits; u = avg + c*sqrt(log(Np)/ch.visits); end
                if u>bestVal, bestVal=u; bestKid=k; end
            end
            curIdx = bestKid;
        end
    
        %% 2‑Expansion
        cur = tree(curIdx);
        if cur.state~=endNode
            neigh = find(~isinf(lambdaMat(cur.state,:)));
            avail = setdiff(neigh, cur.path);
            if isempty(avail)
                tree(curIdx).isFull = true; newIdx = curIdx;
            else
                nxt = avail(randi(numel(avail)));
                new.state=nxt; new.path=[cur.path,nxt];
                new.depth=cur.depth+1; new.visits=0; new.totalLg=0;
                new.parent=curIdx; new.isFull=false;
                tree(end+1)=new; newIdx=numel(tree);
            end
        else
            newIdx = curIdx;
        end
    
        %% 3‑Simulation
        sim.state  = tree(newIdx).state;
        sim.path   = tree(newIdx).path;
        sim.depth  = tree(newIdx).depth;
        lgSum      = 0;
        while sim.state~=endNode
            % (a) direct jump
            if rand<p_direct && ~isinf(lambdaMat(sim.state,endNode))
                lgSum = lgSum + log(p_direct);
                sim.state=endNode; sim.depth=sim.depth+1; sim.path(end+1)=endNode; break; end
            % (b) normal move
            neigh = find(~isinf(lambdaMat(sim.state,:)));
            avail = setdiff(neigh, sim.path);
            if isempty(avail), lgSum = NEG_INF; break; end
            w = lambdaMat(sim.state, avail);  w(isinf(w)) = 0;
            if all(w==0), w = ones(size(w)); end
            p = w / sum(w);
            j = avail(find(cumsum(p) >= rand, 1));
            t  = sim.depth + 1;    lam = lambdaMat(sim.state,j);
            lg = t*log(lam) - lam - gammaln(t+1);   % Poisson log‑pmf (k=1)
            lgSum = lgSum + lg;
            sim.state = j; sim.depth = t; sim.path(end+1) = j;
        end
    
        %% 4‑Backpropagation
        idx = newIdx;
        while idx~=0
            tree(idx).visits   = tree(idx).visits + 1;
            tree(idx).totalLg  = tree(idx).totalLg + lgSum;
            idx = tree(idx).parent;
        end
    
        %% 5‑Update best
        if sim.state==endNode && lgSum>bestLg, bestLg = lgSum; bestPath = sim.path; end
    
        %% 6‑Logging
        log_record(it)  = lgSum;
        prob_record(it) = exp(lgSum);
        step_record(it) = numel(sim.path)-1;
        best_record(it) = bestLg;
        size_record(it) = numel(tree);
    end
    
    %% ========= POST‑PROCESS & SAVE ============
    fprintf('\n最佳路徑: '); disp(bestPath);
    fprintf('log(乘積機率) = %.4g\n', bestLg);
    
    pct = [10 25 50 75 90];
    P   = zeros(numel(pct), iterations);
    for k = 1:iterations
        tmp = log_record(1:k);
        tmp(~isfinite(tmp)) = NEG_INF;      % 轉 finite
        if all(tmp==NEG_INF)
            P(:,k) = NEG_INF;
        else
            P(:,k) = prctile(tmp, pct)';
        end
    end
    x = 1:iterations;
    
    save('mcts2_poisson.mat','lambdaMat','bestPath','bestLg', ...
         'log_record','prob_record','step_record','best_record', ...
         'size_record','P','pct');
    
%% -------- percentile ribbon plot · Steps ---------
    fig = figure('Name','MCTS Progress – Steps','Color','w'); hold on;

    pct  = [10 25 50 75 90];              % 百分位
    Pstp = NaN(numel(pct), iterations);   % 儲存分位數

    for k = 1:iterations
        tmp = step_record(1:k);
        valid = tmp(isfinite(tmp));       % 只保留有限實數
        if ~isempty(valid)
            Pstp(:,k) = prctile(valid, pct)';   % 10/25/…/90 %
        end
    end

    % 左軸：步數帶狀
    fill([x fliplr(x)], [Pstp(1,:) fliplr(Pstp(5,:))], ...
        [0.85 1 0.85], 'EdgeColor','none');          % 10–90 %
    fill([x fliplr(x)], [Pstp(2,:) fliplr(Pstp(4,:))], ...
        [0.55 0.9 0.55], 'EdgeColor','none');        % 25–75 %
    % plot(x, Pstp(3,:), 'k','LineWidth',1.3);          % Median
    % plot(x, step_record, '--g','LineWidth',1.2);      % 每回合實際步數

    legend({'10–90%','25–75%','Median','Steps'}, ...
        'Location','northwest');
    xlabel('Iteration');
    ylabel('Steps');
    grid on;

    saveas(fig, 'mcts2_steps_progress.png');
    end
    