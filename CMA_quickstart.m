% ------------------------------------------------------------------------
%                    http://arxiv.org/abs/1512.01435
% ------------------------------------------------------------------------
% The following script shows how to implement CMA.m with the various
% connection procedures. Written by Martin Ritchie, 2016, University of
% Sussex.
% Dependencies: clustering.m
% ------------------------------------------------------------------------
% Decide on a degree sequence,
D = ones(1,500)*4;
% D = binornd(4,1/125,1,500); 
% D = poissrnd(4,1,500); 

% Creating a random network, i.e., a network using only lines or C2
% subgraphs. Note that (1) the CMA does not require lines as input and (2)
% we now have to specify the sequences for complete squares or 'C4'. 
[Sd, sg] = CMA(D, 'zeros(1,500)','C4');
[A, Sd] = Connect_repeated(Sd, sg);
[C, triangles, triples] = clustering(A)

% Using the same sequence but now with every node incident to a C4 
% subgraph
[Sd, sg] = CMA(D, 'ones(1,500)','C4');
tic
[A, Sd] = Connect_repeated(Sd, sg);
toc
[C, triangles, triples] = clustering(A)

% The number of triangles inputted into the connection procedure, each
% nodes is incident to one C4 corner which is composed of 3 triangles and
% each triangle is counted 3 times, so the number of input triangles is
% given by
Tin = sum(Sd(2,:))
% clustering.m will count each triangle 6 times
Tout = triangles/6

% The number of edges that do not appear in a C4 subgraph is given by
Tin = sum(Sd(1,:))

% We now use the same degree and hyperstub sequences but realised with a
% different algorithm, the unbiased refusal algorithm
tic
[A, Sd, refuse] = Connect_refuse(Sd, sg);
while refuse>0
    [A, Sd, refuse] = Connect_refuse(Sd, sg);
end
toc
% Note that Connect_refuse.m, whilst being unbiased, takes considerably
% longer to run. Its running time will quikly grow as the average degree
% increases. It is the only connection procedure that needs to be run
% within a while loop in order to garuntee that a network will be returned.
D = ones(1,500)*5;
[Sd, sg] = CMA(D, 'ones(1,500)','C4');
tic
[A, Sd, refuse] = Connect_refuse(Sd, sg);
while refuse>0
    [A, Sd, refuse] = Connect_refuse(Sd, sg);
end
toc