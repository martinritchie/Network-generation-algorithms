function A = Connect_erased(Sd, sg)
% ------------------------------------------------------------------------
%                    http://arxiv.org/abs/1512.01435
% ------------------------------------------------------------------------
% The following code will connect the hyperstub sequence(s), Sd, following
% a configuration model type connection procedure. In this version once a
% tuple of nodes is selected the algorithm will form connections regardless
% of self or multi-edges. At the end of the connection process self-edges
% are deleted and multi-edges are collapsed back into a single edge. This
% process has the advantage of  very fast running times but will not
% strictly preserve the degree sequence. Please see the above link for more
% details. Written by Martin Ritchie, University of Sussex, 2016.
% ------------------------------------------------------------------------
% Sd: hyperstub sequences. Sd(1,:) is always lines.
N = length(Sd);
% M: number of subgraphs,
M = length(sg);
% positions(i): the number of hyperstubs in subgraph i,
positions = zeros(1,M);
% ic_subgraph: vector, takes the value 1 if incomplete and 0 otherwise,
ic_subgraph = zeros(1,M);
% corners: each cell entry is a vector of the subgraphs' corner degrees,
corners = cell(1,M);
% edges{i}: the degree sequence of subgraph,
edges = cell(1,M);
% bal{i}: contains the proportions of hyperstubs that comprise subgraph i,
bal = cell(1,M);
% sg{i}: the adjacency matrix of subgraph i,
index = 1;
id = 1;
for i = 1:M
    edges{i} = sum(sg{i});
    corners{i} = unique(edges{i});
    % If the subgraph is not a line and incomplete.
    if mean(edges{i}) >1 && mean(edges{i}) ~= length(edges{i})-1
        ic_subgraph((i)) = 1;
        positions(i) = length(corners{i});
        corner_id(index:index + positions(i)-1) = id;
        bal{i} = hist(edges{i},positions(i));
        index = index + positions(i);
        id = id + 1;
    else
        positions(i) = 1;
        bal{i} = length(edges{i});
        corner_id(index) = id;
        index = index + 1;
        id = id + 1;
    end
    
end
clear index i
%% Hand shake lemma for hyperstub sequences
% Rather than delete surplus hyperstubs that cannot be paired they are
% re-allocated to lines. This preserves the degree sequence,
for i = 2:M
    for j = 1:positions(i)
        % c_pos: used to generate indicies, I, that vary in length
        % depending on the hyperstub types for a given input. For
        % example, in a network composed of lines, triangles and
        % diagonal squares 'I' would need to take the vales: I = 1, 2,
        % 3:4 respectively,
        c_pos = cumsum(positions);
        if positions(i)>1
            if i == 1
                I = 1:c_pos(i);
            else
                I  = c_pos(i-1)+1:c_pos(i);
            end
            % card_comparison: is the value we use to compare the
            % number of hyperstubs,
            card_comparison = sum(Sd(I,:),2);
            % after the following line we require all entries of
            % card_comparison to be equal,
            card_comparison = floor(card_comparison ./ bal{i}');
            [lb, ~] = min(card_comparison);
            [ub, Iub] = max(card_comparison);
            Iub1 = I(Iub);
            while abs(ub-lb) > 0
                % If the numbers are not equal a single surplus
                % hyperstub is selected at random and decomposed back
                % into lines,
                candidate = find(Sd(Iub1,:));
                Ir = randi(length(candidate));
                candidate = candidate(Ir);
                Sd(Iub1,candidate) = Sd(Iub1,candidate) - 1;
                Sd(1,candidate)  = Sd(1,candidate) + corners{i}(Iub);
                card_comparison = sum(Sd(I,:),2);
                card_comparison = floor(card_comparison ./ bal{i}');
                [lb, ~] = min(card_comparison);
                [ub, Iub] = max(card_comparison);
            end
        end
        
    end
end
pcount = 2;
for i = 2:length(positions)
    bcount = 1;
    for j = 1:positions(i)
        % The following ensures the lists are generate with correct
        % multiplicty. Do not check lines, i~=1, until after we have
        % decomposed excess corners back into lines
        while mod(sum(Sd(pcount,:)),bal{i}(bcount))>0
            index = find(Sd(pcount,:),1);
            Sd(pcount,index) = Sd(pcount,index) - 1;
            Sd(1,index) = Sd(1,index) + corners{i}(j);
        end
        pcount = pcount + 1;
        bcount = bcount + 1;
    end
end
if mod(sum(Sd(1,:)),2)~=0
    index = find(Sd(1,:),1);
    Sd(1,index) = Sd(1,index) - 1;
end
%% Forming dynamic lists for the HCM process
Sg_bin = cell(1,sum(positions));
index1 = 1;
for i = 1:M
    for j = 1:positions(i)
        index2 = 1;
        for k = 1:N
            for l = 1:Sd(index1,k)
                Sg_bin{index1}(index2) = k;
                index2 = index2 + 1;
            end
        end
        index1 = index1 + 1;
    end
end
A = zeros(N);
%% Connection process
for i = 1:M
    c = sum(sg{i});
    % Generate the corect proportions of corners.
    if ic_subgraph(i) == 1;
        porportions =  hist(c, corners{i});
        porportions = porportions(porportions~=0);
    else
        porportions = length(sg{i});
    end
    % Check for non-empty bins. Sg_bin is a vector of cells. Each
    % entry is a dynamic list of nodes.
    while ~isempty(Sg_bin) && sum(cellfun(@isempty,Sg_bin(1:positions(i))))==0
        % (Initialise) nodei, candidate node indices from bin.
        nodei = cell(1,positions(i));
        % (Initialise) nodes, the nodes labels.
        nodes = nodei;
        % Initial selection:
        for j = 1:positions(i)
            nodei{j} = randi(length(Sg_bin{j}),1,porportions(j));
            nodes{j} = Sg_bin{j}(nodei{j});
        end
        % This line connects the subgraph
        A([nodes{:}],[nodes{:}]) = sg{i};
        % Then remove the entries from the bin
        for j = 1:positions(i)
            Sg_bin{j}(nodei{j}) = [];
        end
        test  = sum(cellfun(@isempty,Sg_bin(1:positions(i))));
        if test > 0
            for j = 1:positions(i)
                Sg_bin(1) = [];
            end
            break
        end
    end
    if  isempty([Sg_bin{:}])
        for ii = 1:N
            A(ii,ii)=0;
        end
        A = sparse(A);
        return
    end
end