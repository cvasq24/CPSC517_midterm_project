function [P] = ap_min_degree_nodense(Amtx)
    n = size(Amtx,1);
    
    % We keep track of the following:
    % - head: head of degree bucket linked list. E.g., head(1) contains
    % nodes with degree 0, head(2) nodes w/ degree 1, and so on.
    % - next: pointer
    % - last: previous pointer
    head = zeros(n+1, 1);
    next = zeros(n, 1);
    last = zeros(n, 1);

    % - nv: supervariable sizes
    % - d: approximate degrees   
    % - w: weights for Algorithm 2
    % - wflg: flag to avoid clearing w
    nv = ones(n, 1);
    d = zeros(n, 1);
    w = zeros(n, 1);
    wflg = 2; % first iteration, nodes have starting weight of 1
    
    % - A: variable's adjacent variables
    % - E: variable's adjacent elements
    % - Le: element's adjacent variables
    % - sv_members: member of the supervariable with principal variable i
    A = cell(n, 1);
    E = cell(n, 1);
    Le = cell(n, 1);
    sv_members = cell(n, 1);
    
    % Initialize: A[i] are a node's adjacent variables, d[i] is |Ai|.
    % Each node is its own supervariable, hence sv_members[i] = i.
    for i = 1:n
        A{i} = setdiff(find(Amtx(i,:)), i);
        d(i) = length(A{i});
        sv_members{i} = i;
    end
    
    % Dense row detection
    % Removes nodes (rows) with more non-zeros (i.e. higher degree) than
    % threshold (i.e. set the size of its supervariable to 0).
    dense_rows = [];

    % Populate the initial degree buckets (by degree value).
    % For the sparse rows, detach dense nodes from A[i], update degree
    % d[i], and insert into the linked list.     
    for i = 1:n
        if nv(i) > 0 
            A{i} = setdiff(A{i}, dense_rows);
            d(i) = length(A{i});
            
            % Insert node i into its degree bucket
            % Get its degree. If there is another node at the front of the
            % degree bucket (inext), then link them so i is in front of
            % inext. 
            % Finally, make node i the new head of the bucket
            deg = d(i) + 1; % 1-indexed
            inext = head(deg);
            if inext ~= 0, last(inext) = i; end
            next(i) = inext;
            head(deg) = i;
            last(i) = 0;
        end
    end

    k = 1; 
    P = []; 
    mindeg = 1; % skip bucket with degree 0
    total_pivots = n - length(dense_rows);
    fprintf('Crunching matrix... Progress: 000%%');
    while k <= n - length(dense_rows)
        % 1. Select pivot. O(1) b/c of the bucket approach
        % Iterate over degree buckets 1, 2, ... until the bucket is
        % non-empty, grab its head - that is our pivot
        p = 0;
        for deg = mindeg:(n+1)
            p = head(deg);
            if p ~= 0, break; end
        end
        if p == 0, break; end
        mindeg = deg; % update minimum degree for next iteration

        % 2. Remove pivot from the bucket & update arrays
        % next(p) is the next node w/ same degree as p, so grab it as next
        % next head of bucket w/ head(deg)
        % (e.g. if nodes (2,4) have degree=3, head(3) = 2, but
        % after we set p=2, head(3) = 4. 
        inext = next(p);
        if inext ~= 0, last(inext) = 0; end
        head(deg) = inext;
        nv_piv = nv(p); 
        nv(p) = 0; % mark as 'dead' node (not a variable anymore)

        % 3. Construct Lp (equation 3.1 in AMD paper)
        % Get its adjacent variable, and the neighbours of its adjacent
        % elements. 
        Lp = A{p};
        Ep_arr = reshape(E{p}, 1, []);
        if ~isempty(Ep_arr)
            % MATLAB Magic: [Le{Ep_arr}] automatically unwraps and concatenates 
            % all the arrays instantly without a for-loop.
            Lp = [Lp, Le{Ep_arr}]; 
        end
        Lp = unique(Lp);
        Lp = Lp(nv(Lp) > 0); 
        Lp(Lp == p) = []; 
        Le{p} = Lp;

        % 4. Algorithm 2
        for e = reshape(E{p}, 1, [])
            w(e) = 0;
        end

        if mod(k, 1000) == 0
            fprintf('\b\b\b\b%03d%%', round((k / total_pivots) * 100));
        end
        
        % Scan the elements for i in Lp
        % Iterate over p's adjacent variables' element neighbours (e)
        % This is 'weight' w(e). 
        % We use the wflg here to help speed this step up. Taken from the
        % Fortran implementation. This avoids re-setting w(e) = zeros(n,1)
        % every time we select a new pivot. 
        for i = reshape(Lp, 1, [])
            for e = reshape(E{i}, 1, [])
                if e == p, continue; end
                if w(e) < wflg
                    % Get 'non-dead' variables adjacent to e
                    % Weight is as in algorithm 2: |Le|
                    alive_Le = Le{e}(nv(Le{e}) > 0);
                    w(e) = sum(nv(alive_Le)) + wflg; 
                end
                % As in algorithm 2: w(e) = w(e) - |i|
                w(e) = w(e) - nv(i);
            end
        end

        % 5. Updating degrees & absorptions into pivot
        for i = reshape(Lp, 1, [])
            if nv(i) == 0, continue; end 
            
            % The pivot p is 'removed' from graph, so variable i has to be 
            % updated in linked list. This matters for d(i) computation.
            % Later on variable i is re-inserted in the right bucket.

            % Unhook i from its old degree bucket and link the nodes before
            % and after it in the linked list (i.e. prev->nxt). 
            old_deg = d(i) + 1;
            prev = last(i); 
            nxt = next(i);
            if prev ~= 0, next(prev) = nxt; else, head(old_deg) = nxt; end
            if nxt ~= 0, last(nxt) = prev; end

            % Element absorption as in Algorithm 1
            % Also grab the third part in the third term of the min
            % operation in equation 4.1 (approx degree) in one go
            % Unlike 4.1, where algorithm checks if w(e) >= 0, here we
            % dont worry about that because previously we scanned all
            % elements in Ei 

            Ei_arr = reshape(E{i}, 1, []);
            Ei_arr(Ei_arr == p) = []; % Remove p instantly
            
            % Create a logical mask of which elements survive
            valid_mask = w(Ei_arr) > wflg; 
            new_Ei = Ei_arr(valid_mask);
            
            % Calculate sum_we in one shot
            sum_we = sum(w(new_Ei) - wflg);
            
            E{i} = [new_Ei, p];
            
            % Update adjacency: remove redundant entires 
            Ai = A{i}(nv(A{i}) > 0); % only 'non-dead' variables
            Ai = setdiff(Ai, Lp); % as in Algorithm 1
            A{i} = Ai;

            % Approx Degree (Equation 4.1)
            term_Ai = sum(nv(Ai)); 
            term_Lp = sum(nv(Lp)) - nv(i);
            d_new = term_Ai + term_Lp + sum_we;
            
            % Mass Elimination
            if d_new == 0
                nv_piv = nv_piv + nv(i);
                sv_members{p} = [sv_members{p}, sv_members{i}];
                nv(i) = 0; % Dead
            else
                d(i) = min(n - k, d_new);
                % Re-hook i into new degree bucket
                new_deg = d(i) + 1;
                inext = head(new_deg);
                if inext ~= 0, last(inext) = i; end
                next(i) = inext; head(new_deg) = i; last(i) = 0;
                if new_deg < mindeg, mindeg = new_deg; end
            end
        end

        % 6. Supervariable detection via hash function
        Lp_alive = Lp(nv(Lp) > 0);
        if length(Lp_alive) > 1
            h = zeros(length(Lp_alive), 1);
            for idx = 1:length(Lp_alive)
                v = Lp_alive(idx);
                % Use hash function from AMD paper (page 897)
                h(idx) = mod(sum(A{v}) + sum(E{v}), n-1) + 1; 
            end
            
            [unique_h, ~, ic] = unique(h);
            for u_idx = 1:length(unique_h)
                nodes = Lp_alive(ic == u_idx);
                if length(nodes) > 1
                    for idx1 = 1:length(nodes)
                        i = nodes(idx1); if nv(i) == 0, continue; end
                        for idx2 = idx1+1:length(nodes)
                            j = nodes(idx2); if nv(j) == 0, continue; end
                            
                            if isequal(sort(A{i}), sort(A{j})) && isequal(sort(E{i}), sort(E{j}))
                                nv(i) = nv(i) + nv(j);
                                sv_members{i} = [sv_members{i}, sv_members{j}];
                                nv(j) = 0; % Mark j as dead
                                
                                % Unhook j from bucket
                                old_deg = d(j) + 1;
                                prv = last(j); nxt = next(j);
                                if prv ~= 0, next(prv) = nxt; else, head(old_deg) = nxt; end
                                if nxt ~= 0, last(nxt) = prv; end
                                
                                A{j} = []; E{j} = [];
                            end
                        end
                    end
                end
            end
        end

        % 7. Finalize pivot: clear its A/E arrays, insert it into P, update
        % k = k + |p| as in Algorithm 1
        P = [P, sv_members{p}]; 
        k = k + nv_piv; 
        A{p} = []; 
        E{p} = [];
        
        % Increment wflg so we never have to zero out the w array
        wflg = wflg + n; 
    end  
    
    % Insert dense rows that we excluded before as last
    P = [P, dense_rows];
    
    % Post-ordering
   [~, ~, ~, post] = symbfact(Amtx(P, P));
    P = P(post);
end