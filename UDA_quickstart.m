% ------------------------------------------------------------------------
%                    http://arxiv.org/abs/1512.01435
% ------------------------------------------------------------------------
% The following script shows how to implement UDA.m with the various
% connection procedures. Written by Martin Ritchie, 2016, University of
% Sussex.
% Dependencies: clustering.m
% ------------------------------------------------------------------------
% Decide on a degree sequence,
D = ones(1,500)*4;
% D = binornd(5,1/100,1,500); 
% D = poissrnd(5,1,500); 

% Creating a random network, i.e., a network using only lines or C2
% subgraphs.
[Sd, sg] = UDA(D, 'C2');
[A, Sd] = Connect_repeated(Sd, sg);
[C, triangles, triples] = clustering(A)

% Using the same sequence but now with triangle, or C3, subgraphs
[Sd, sg] = UDA(D, 'C2','C3');
tic
[A, Sd] = Connect_repeated(Sd, sg);
toc
[C, triangles, triples] = clustering(A)

% The number of triangles inputted into the connection procedure, diving by
% three to eliminate the triple count
Tin = sum(Sd(2,:))/3
% clustering.m will count each triangle 6 times
Tout = triangles/6

% We now use the same degree and hyperstub sequences  but realised with a
% different algorithm, the unbiased refusal algorithm
tic
[A, Sd, refuse] = Connect_refuse(Sd, sg);
while refuse>0
    [A, Sd, refuse] = Connect_refuse(Sd, sg);
end
toc
% Note that Connect_refuse.m, whilst being unbiased, takes considerably
% longer to run. Its running time will quikly grow as the average degree
% increases.
D = ones(1,500)*5;
[Sd, sg] = UDA(D, 'C2','C3');
tic
[A, Sd, refuse] = Connect_refuse(Sd, sg);
while refuse>0
    [A, Sd, refuse] = Connect_refuse(Sd, sg);
end
toc

