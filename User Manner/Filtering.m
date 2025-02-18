function[E]=Filtering(Cluster,W)
L=length(Cluster);
E(length(W))=0;
for i=1:L
    Condidate=Cluster(i).cluster;
    W_cluster(i)=mean(W(Condidate));
end
k=find(W_cluster>0.65*max(W));
[a,b]=sort(W_cluster(k));
if length(b)>0
for j=1:length(b)
    A=Cluster(k(b(j))).cluster;
    [a1,b1]=sort(W(A));
    E(A(b1(length(b1))))=1;
    if length(b1)>1
        E(A(b1(length(b1)-1)))=1;
    end
    if length(b1)>2
        E(A(b1(length(b1)-2)))=1;
    end 
end
end
    