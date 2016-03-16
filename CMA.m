function [Sd, sg] = CMA(D, varargin)
% ------------------------------------------------------------------------
%                    http://arxiv.org/abs/1512.01435
% ------------------------------------------------------------------------
% The following code will generate a complex network following the
% cardinality matching algorithm (CMA). Please see the above link for more
% details. Written by Martin Ritchie, University of Sussex, 2016.
%% Inputs / output
% D: Degree sequence,
% varargin: a list of subgraph sequences followed by the list of subgraphs,
% A: adjacency matrix (sparse format).
% Sd: hyperstub sequence(s)
%% Dependencies
% * subgraphs.mat
%% Example call
%  D = ones(1,1000)*6;
%  A = CM_card(D, 'ones(1,1000)','ones(1,1000)', 'C3','Toast');
%% Complete/incomplete subgraph classification
% This only needs to be performed once for each subgraphs specified in
% varargin. Complete subgraphs do not require any identification, all nodes
% are equivalent. Incomplete subgraphs require that each and every node be
% given a unique label and each label has its stub cardinality associated
% with it. This section is common to the code for both the UDA and CMA.
% N: number of nodes.
N = length(D);
% NN: a node list used for logical indexing.
NN = 1:N;
% M: number of subgraphs.
M = floor(length(varargin)/2);
load('subgraphs',varargin{M+1:end});
% positions(i): the number of hyperstubs in subgraph i.
positions = zeros(1,M);
% ic_subgraph: vector, takes the value 1 if incomplete and 0 otherwise.
ic_subgraph = zeros(1,M);
% corners: each cell entry is a vector of the subgraphs' corner degrees.
corners = cell(1,M);
% edges{i}: the degree sequence of subgraph.
edges = cell(1,M);
% bal{i}: contains the proportions of hyperstubs that comprise subgraph i.
bal = cell(1,M);
% sg{i}: the adjacency matrix of subgraph i.
sg = cell(1,M);
index = 1;
id = 1;
for i = 1:M
    sg{i} = eval(varargin{M+i});
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
clear index
% Ensure that lines are included.
if  ~isequal(ones(2) - eye(2),sg{1})
    sgt = sg;
    sg{1} = ones(2) - eye(2);
    sg(2:length(sgt)+1) = sgt;
    clear sgt
end
%% The CMA process
% The following forms, 'S_card_v', a single vector resulting from
% concatenating all hyperstub sequences. Each each hyperstub
% represented in this list by its stub cardinality multiplied by the
% edge count of the hyperstub. card_error: returns '0' if it can pair
% all hyperstub degrees to nodes without any problems.
%% Converting the input subgraph sequences into hyperstub sequences
% S: subgraph degree sequence(s).
S = cell(1,M);
% S_card: A sequence specifying the total number of edges required
% by a given subgraph degree.
S_card = S;
for i = 1:M
    S{i} = eval(varargin{i});
    S_temp = S;
    if ic_subgraph(i)==1
        % p: the proportions that a hyperstub appears within a
        % subgraph.
        p = hist(edges{i},length(unique(corners{i})));
        p = p/sum(p);
        for k = 1:N
            % For incomplete subgraphs we multinomially divide the
            % subgraph degree between the hyperstubs.
            S{i}(1:positions(i),k) = mnrnd(S_temp{i}(k),p)';
            % S_card: stores the required stub cardinality induced
            % by the hyperstub degree.
            S_card{i}(1:positions(i),k) = corners{i}'.*S{i}(1:positions(i),k);
        end
        test = sum(S{i}(1:positions(i),:),2);
        while length(unique(test))~=1
             S_temp{i} = eval(varargin{i});
            for k = 1:N
                S{i}(1:positions(i),k) = mnrnd(S_temp{i}(k),p)';
                S_card{i}(1:positions(i),k) = corners{i}'.*S{i}(1:positions(i),k);    
            end
            test = sum(S{i}(1:positions(i),:),2);
        end
    else
        % For complete subgraph.
        S_card{i} = S_temp{i}.*corners{i};
    end
end
% The following re-compiles all of S_card and S_card_ic into a
% vector. This vector is compared to the degree sequence to find
% suitable pairings between nodes and hyperstub degrees.
SM = zeros(sum(positions),N);
ind = 1;
for i = 1:M
    [d1, ~] = size(S_card{i});
    SM(ind:ind+d1-1,:) = S_card{i};
    ind = ind + d1;
end
d = D;
SMt = sum(SM,1);
[SMt, I] = sort(SMt);
SM = SM(:,I);
Sd =  zeros(sum(positions),N);
n = N;
while n > 0
    eligible = NN(d>=SMt(n));
    d1 = length(SM(:,n));
    % Since the hyperstub sequences are not ordered the following
    % reorders them, a single hyperstub at a time, until the degree
    % induced by the hypersyubs will fit into what is remaining of the
    % classical degree sequence.
    while isempty(eligible)
        corner_card = zeros(1,M);
        for i = 1:d1
            corner_card(corner_id(i)) =  corner_card(corner_id(i)) + SM(i,n);
        end
        % Prune a hyperstub at a time. Taking the smallest possible
        % reduction each time.
        [~, min_index_1] = max(corner_card);
        min_index_2 = find(corner_id==min_index_1);
        
        temp = double(SM(min_index_2,n)>0);
        temp = temp.*corners{min_index_1}';
        SM(min_index_2,n) = SM(min_index_2,n) - temp;
        
        % Now reallocate to another node.
        new_home = NN(SMt(1:(n-1)) < sum(SM(:,n)));
        new_home = new_home(randi(length(new_home)));
        while  sum(SM(min_index_2,new_home)) > 0
            new_home = NN(SMt(1:(n-1)) < sum(SM(:,n)));
            new_home = new_home(randi(length(new_home)));
        end
        SM(min_index_2,new_home) = SM(min_index_2,new_home) + temp;
        SMt(n) = sum(SM(:,n));
        SMt(new_home) = sum(SM(:,new_home));
        
        eligible = NN(d>=SMt(n));
        
        %------------------------------------------------------
    end
    candidate = eligible(randi(length(eligible)));
    Sd(:,candidate) = Sd(:,candidate) + SM(:,n)./[corners{:}]';
    % The same node may not be allocated the same type of corner more than
    % once. 
    d(candidate) = 0;
    n = n - 1;
end
V1 = repmat([corners{:}]',1,N);
d = D - sum(Sd.*V1,1);
Sd = [d; Sd];
end