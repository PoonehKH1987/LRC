function[Cluster]=MCL(W_graph,Vertices,n)
%Adj=W_graph.Adj;
%V=W_graph.Vertices;
[N]=Normalization(W_graph);
for i=1:n
    M=N*N;
    [N]=Inflation2(M);
end

[Cluster]=MCLCluster(N,Vertices);