function [Sd, sg] = UDA(D, varargin)
% ------------------------------------------------------------------------
%                    http://arxiv.org/abs/1512.01435
% ------------------------------------------------------------------------
% The following code will generate a complex network following the
% Underdetermined Diophantine Algorithm (UDA). Please see the above link
% for more details. Written by Martin Ritchie and Istvan Kiss, University
% of Sussex, 2016.
%% Inputs / outputs
% D: degree sequence,
% varargin: subgraph names (string),
% Sd: subgraph sequences,
% Sg: cell array of subgraph adjacency matrices. 
%% Dependencies
% * dio_recur.m
% * subgraphs.mat
%% Example call
% [Sd, sg] = UDA(D, 'C2','C3','Toast');
%% Complete/incomplete subgraph classification
% This only needs to be performed once for each subgraphs specified in
% varargin. Complete subgraphs do not require any identification, all nodes
% are equivalent. Incomplete subgraphs require that each and every node be
% given a unique label and each label has its stub cardinality associated
% with it.
N = length(D);
load('subgraphs',varargin{:});
M = length(varargin);
% ic_subgraph: incomplete subgraph vector, takes the value 1 if incomplete
% and 0 otherwise.
ic_subgraph = zeros(1,M);
% corners: each cell entry is a vector of the degree sequences for the
% subgraphs.
corners = cell(1,M);
% positions: the number of nodal positions required for the given set of
% subgraphs.
positions = zeros(1,M);
% Balance{i}(:) gives a vector of proportions for each corner type.
bal = cell(1,M);
edges = cell(1,M);
% sg{i}: the adjacency matrix of subgraph i.
sg = cell(1,M);
for i = 1:M
    sg{i} = eval(varargin{i});
    edges{i} = sum(eval(varargin{i}));
    if mean(edges{i}) >1 && mean(edges{i}) ~= length(edges{i})-1
        ic_subgraph(i) = 1;
        corners{i} = unique(sum(eval(varargin{i})));
        % For incomplete subgraphs composed of g nodes, we use g positions.
        positions(i) = length(corners{i});
        bal{i} = hist(edges{i},length(unique(corners{i})));
    else
        % For complete subgraphs we use a single position:
        positions(i) = 1;
        corners{i} = unique(sum(eval(varargin{i})));
        bal{i} = length(edges{i});
    end
end
%% Creating the Diophantines' solution space(s)
% solution_space: contains all solutions to the under-determined equations
% that are used to establish stub configurations. solution_space{i}
% corresponds to the solution space of a node degree, k=i.
solution_space = cell(1,max(D));
% classes: each class, classes(i), represents a solution space for
% degree i.
classes = zeros(1,max(D));
solution_class = zeros(1,length(max(D)));
k = 1;
for i = 1:max(D)
    solution =  dio_recur([corners{:}],i);
    solution_space{i} =  solution';
    clearvars -global sol_space
    [~, x2] = size(solution_space{i});
    classes(i) = i;
    for j = 1:x2
        solution_class(k) =i;
        k = k + 1;
    end
end
%% Creating degree sequences
Sd = zeros(sum(positions),N);
for i = 1:N
    if D(i)>0
        [~, l] = size(solution_space{D(i)});
        r = randi(l);
        Sd(:,i) = solution_space{D(i)}(:,r);
    end
end
end