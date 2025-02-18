function[N]=Inflation2(M)
%M=M.^2;
[m,n]=size(M);
N(m,n)=0;
for i=1:n
    A=(M(:,i).*M(:,i))./(sum(M(:,i).*M(:,i)));
    N(:,i)=A;
end
    