function[W,W_graph]=weigthed_Gragh(Un_W_graph,X)
W=LR(X);
L=length(W);

I(L,L)=0;
for i=1:L
    I(i,i)=1;
end


Un_W_graph=Un_W_graph+I;

for i=1:L
    Neghboor(i).N=find(Un_W_graph(i,:)>0);
end

Wgraph(L,L)=0;
for i=1:L
        N_i=Neghboor(i).N;
    for j=i+1:L
        N_j=Neghboor(j).N;
        N_ij=intersect(N_i,N_j);
        Wgraph(i,j)=(sum(W(N_ij)))/(sum(W(N_i))+sum(W(N_j))-sum(W(N_ij)));
        Wgraph(j,i)=Wgraph(i,j);
    end
end
W_graph=Wgraph;
