function [P] = min_degree(Amtx)
    n = size(Amtx,1);
    
    A = cell(n,1);
    E = cell(n,1);
    d = zeros(n,1);
    sv = cell(n,1); % super variables

    % ---- Initialize A/E/d/i  ----
    for i = 1:n
        A{i} = setdiff(find(Amtx(i,:)), i); % exclude itself
        E{i} = [];
        d(i) = length(A{i});
        sv{i} = i; % initialized as themselves
    end 

    k = 1;    
    P = [];
    V = 1:n;
    Vbar = [];
    Le =  cell(n,1); % for each element
    
    while (k <= n)
        % ---- Mass eliminiation ----
        % 1. Find current pivot (that minimizes degree) from variables
        [~, min_idx] = min(d(V));
        p = V(min_idx);
        
        % 2. Get its Lp using 3.1 from page 3 in AMD paper
        Lp = A{p};
        for e = E{p}
            Lp = union(Lp, Le{e});
        end
        Lp = setdiff(Lp, sv{p});
        Le{p} = Lp;
        
        for i = Lp
            A{i} = setdiff(A{i}, [Lp, p]);
            E{i} = union(setdiff(E{i}, E{p}), p);
            d(i) = get_external_degree(A,E,sv,Le,i); % external degree
        end

        % ---- Supervariable detection ----
        % Iterate (brute force) over (i,j) supervariable pairs in Lp
        % If indistinguishable, remove j
        removed = false(n, 1);
        for ii = 1:length(Lp)
            i = Lp(ii);
            if removed(i), continue; end

            for jj = ii+1:length(Lp)
                j = Lp(jj);
                if removed(j), continue; end

                if is_pair_indist(A,E,i,j) == true
                    sv{i} = union(sv{i}, sv{j});
                    V = setdiff(V, j);
                    d(i) = d(i) - length(sv{j});
                    A{j} = [];
                    E{j} = [];
                    removed(j) = true;
                end
            end 
        end
        
        % ---- Convert pivot onto an element ----
        Vbar = setdiff(union(Vbar, p), E{p});
        V = setdiff(V, p);
        A{p} = [];
        E{p} = [];
        k = k + length(sv{p});
        P = [P, sv{p}];
    end  
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

function  bool = is_pair_indist(A,E,i,j)
    same_var_adj = isequal(sort(A{i}), sort(A{j}));
    same_ele_adj = isequal(sort(E{i}), sort(E{j}));
    bool = same_var_adj && same_ele_adj;
end