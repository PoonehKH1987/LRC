function[Protein]=Antigen_Surface_Graph2(Sur_Res)

P=pdbread('3NH7.pdb');
L=length(P.Model.Atom);
k=1;
s=1;
for i=1:L
    if P.Model.Atom(i).chainID==char('A')
        if intersect(P.Model.Atom(i).resSeq,Sur_Res)>0
            S(k)=P.Model.Atom(i).resSeq;
            C(k,1)=P.Model.Atom(i).X;
            C(k,2)=P.Model.Atom(i).Y;
            C(k,3)=P.Model.Atom(i).Z;
            k=k+1;
            if P.Model.Atom(i).resSeq==Sur_Res(1)
                if s==1
                   Res_Name{s}=char(P.Model.Atom(i).resName);
                   s=s+1;
                end
            elseif P.Model.Atom(i).resSeq>P.Model.Atom(i-1).resSeq
                Res_Name{s}=char(P.Model.Atom(i).resName);
                s=s+1;
            end
        end
    end
end
x=C(:,1)';
y=C(:,2)';
z=C(:,3)';

T = delaunay3(x,y,z);
L1=length(T);
Adj=[];
Adj(max(Sur_Res),max(Sur_Res))=0;
i=1;
while i<L1+1
    if Adj(S(T(i,1)),S(T(i,2)))==0
        d=sqrt(sum((C(T(i,1),:)-C(T(i,2),:)).^2));
            if d<=6
               Adj(S(T(i,1)),S(T(i,2)))=1;
               Adj(S(T(i,2)),S(T(i,1)))=1;
            end
    elseif Adj(S(T(i,1)),S(T(i,3)))==0
        d=sqrt(sum((C(T(i,1),:)-C(T(i,3),:)).^2));
            if d<=6
               Adj(S(T(i,1)),S(T(i,3)))=1;
               Adj(S(T(i,3)),S(T(i,1)))=1;
            end
    elseif Adj(S(T(i,1)),S(T(i,4)))==0
        d=sqrt(sum((C(T(i,1),:)-C(T(i,4),:)).^2));
            if d<=6
               Adj(S(T(i,1)),S(T(i,4)))=1;
               Adj(S(T(i,4)),S(T(i,1)))=1;
            end
    elseif Adj(S(T(i,2)),S(T(i,3)))==0
        d=sqrt(sum((C(T(i,2),:)-C(T(i,3),:)).^2));
            if d<=6
               Adj(S(T(i,2)),S(T(i,3)))=1;
               Adj(S(T(i,3)),S(T(i,2)))=1;
            end
    elseif Adj(S(T(i,2)),S(T(i,4)))==0
        d=sqrt(sum((C(T(i,2),:)-C(T(i,4),:)).^2));
            if d<=6
               Adj(S(T(i,2)),S(T(i,4)))=1;
               Adj(S(T(i,4)),S(T(i,2)))=1;
            end
    elseif Adj(S(T(i,3)),S(T(i,4)))==0
        d=sqrt(sum((C(T(i,3),:)-C(T(i,4),:)).^2));
            if d<=6
               Adj(S(T(i,3)),S(T(i,4)))=1;
               Adj(S(T(i,4)),S(T(i,3)))=1;
            end
    end
    i=i+1;
end

Protein.SurfaceResidue.Number=Sur_Res;
Protein.SurfaceResidue.Name=Res_Name;
Protein.SurfaceResidue.Adjancency=Adj(Sur_Res,Sur_Res);


