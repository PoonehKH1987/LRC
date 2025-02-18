function[Cluster]=MCLCluster(N,Vertices) 

n=length(N);
Cluster=[];
k=1;
for i=1:n
X=find(N(i,:)>0);
x=length(X);
if x>0
Cluster(k).cluster=X;
Cluster(k).Residues=Vertices(X);
k=k+1;
end
end
