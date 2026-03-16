function [P] = ap_min_degree2(Amtx)
    n = size(Amtx,1);
    %  P_exact = min_degree(Amtx);  % run once before the loop
    
    V = 1:n;
    Vbar = [];
    
    A = cell(n,1);          % adjacent variables
    E = cell(n,1);          % adjacent elements
    d = zeros(n,1);         % approximate degree
    sv = cell(n,1);         % supervariables
    sv_size = ones(n,1);    % supervariable sizes
    
    P = [];                 % pivot ordering
    w = -ones(n, 1);        % weights for Algorithm 2
    Le = cell(n,1);         % element variable lists

    % ---- Initialize ----
    for i = 1:n
        A{i} = setdiff(find(Amtx(i,:)), i);
        E{i} = [];
        d(i) = length(A{i});
        sv{i} = i;
    end 

    % Dense row detection (MATLAB default: threshold = 10*sqrt(n))
    dense_threshold = 10* sqrt(n);
    dense_rows = [];
    for i = 1:n
        if d(i) > dense_threshold && d(i) > 16
            dense_rows = [dense_rows, i];
        end
    end
    
    if ~isempty(dense_rows)
        V = setdiff(V, dense_rows);
        for i = V
            A{i} = setdiff(A{i}, dense_rows);
            d(i) = length(A{i});
        end
        % Clear dense rows' adjacency
        for i = dense_rows
            A{i} = [];
            E{i} = [];
        end
    end

    k = 1;    
    while (k <= n)
        % 1. Find pivot with minimum degree
        % [~, min_idx] = min(d(V));
        % p = V(min_idx);
        d_min = min(d(V));
        candidates = V(d(V) == d_min);
        p = candidates(randi(length(candidates)));
        
        % % Then inside the while loop, after picking p:
        % if p ~= P_exact(k) && k <= 1000
        %     fprintf('k=%d: approx picks %d (d=%d), exact would pick %d (d=%d)\n', ...
        %         k, p, d(p), P_exact(k), d(P_exact(k)));
        % end

        % 2. Compute Lp (reachable set of p)
        Lp = A{p};
        for e = E{p}
            Lp = union(Lp, Le{e});
        end
        Lp = setdiff(Lp, sv{p});
        Lp = intersect(Lp, V);  % safety: only active variables
        Le{p} = Lp;

        % 3. Algorithm 2 - compute w(e) = |Le \ Lp| for all elements
        for i = Lp
            for e = E{i}
                if w(e) < 0
                    w(e) = length(Le{e});
                end
                w(e) = w(e) - 1;
            end
        end

        % 4. Aggressive absorption: if Le(e) ⊆ Lp ∪ sv{p}, absorb e
        for e = Vbar
            if w(e) == 0
                for i = Le{e}
                    E{i} = setdiff(E{i}, e);
                end
                Vbar = setdiff(Vbar, e);
                Le{e} = [];
            end
        end

        % 5. Update adjacency and degree for all variables in Lp
        for i = Lp
            A{i} = setdiff(A{i}, [Lp, sv{p}]);
            E{i} = union(setdiff(E{i}, E{p}), p);

            d_old = d(i);
            d(i) = get_approx_degree(A, E, Le, sv_size, w, n, k, i, p, d_old);
        end

        % 6. Supervariable detection with hashing
        h = zeros(n, 1);
        for ii = 1:length(Lp)
            i = Lp(ii);
            h(i) = mod(sum(A{i}) + sum(E{i}), n - 1) + 1;
        end

        removed = false(n, 1);
        for ii = 1:length(Lp)
            i = Lp(ii);
            if removed(i), continue; end

            for jj = ii+1:length(Lp)
                j = Lp(jj);
                if removed(j), continue; end
                if h(i) ~= h(j), continue; end  % hash filter

                if is_pair_indist(A, E, i, j)
                    sv{i} = union(sv{i}, sv{j});
                    sv_size(i) = sv_size(i) + sv_size(j);  % keep sv_size in sync
                    V = setdiff(V, j);
                    d(i) = d(i) - sv_size(j);
                    A{j} = [];
                    E{j} = [];
                    removed(j) = true;
                end
            end 
        end
    
        % 7. Reset weights
        for i = Lp
            for e = E{i}, w(e) = -1; end
        end
        
        % 8. Convert pivot to element
        Vbar = setdiff(union(Vbar, p), E{p});
        V = setdiff(V, p);
        A{p} = [];
        E{p} = [];
        k = k + length(sv{p});
        P = [P, sv{p}];
    end  
    P = [P, dense_rows];
end

function di = get_approx_degree(A, E, Le, sv_size, w, n, k, i, p, d_old)
    % Term 1: trivial upper bound
    term1 = n - k - 1;

    % Term 2: previous degree + new neighbors from Lp
    lp_minus_i = length(Le{p}) - sv_size(i);
    term2 = d_old + lp_minus_i - 1;

    % Term 3: |Ai| + |Lp \ i| + sum of w(e) for other elements
    %   A{i} already excludes sv{i} (self-loops removed at init)
    %   and has been pruned of Lp and sv{p}, so length(A{i}) is clean
    term3 = length(A{i}) + lp_minus_i;

    for e = setdiff(E{i}, p)
        if w(e) >= 0
            term3 = term3 + w(e);
        else
            % Element not scanned by Alg2 — use full |Le| as fallback
            term3 = term3 + length(Le{e});
        end
    end

    di = min([term1, term2, term3]);
end

function di = get_external_degree(A,E,sv,Le,i)
    % First term
    di = length(setdiff(A{i}, sv{i}));

    % Second term
    second_term = [];
    for e = E{i}
        second_term = union(second_term, Le{e});
    end 
    di = di + length(setdiff(second_term, sv{i}));   
end

function bool = is_pair_indist(A, E, i, j)
    same_var_adj = isequal(sort(A{i}), sort(A{j}));
    same_ele_adj = isequal(sort(E{i}), sort(E{j}));
    bool = same_var_adj && same_ele_adj;
end
