function [phi, Triangles, Triples]=clustering(G)
% Computes the global clustering coefficient of G. 
% [phi, Triangles, Triples]=clustering(G)
G2 = G*G; 
G3 = G2*G;


Triangles=full(sum(diag(G3)));

Triples=full(sum(sum(G2))-sum(diag(G2)));

phi = Triangles/(Triples);

end