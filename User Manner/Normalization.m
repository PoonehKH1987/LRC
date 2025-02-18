function[N]=Normalization(W_graph)
n=length(W_graph);
I(n,n)=0;
for i=1:n
    I(i,i)=1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
W_graph=W_graph+I;
S=sum(W_graph);
N(n,n)=0;
for i=1:n
    a=W_graph(:,i)/S(i);
    N(:,i)=a;
end
