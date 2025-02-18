function[Sur_Res,Sur_Atom]=surface_residues(ASA)
A=Importdata(char(ASA));
X=find(A.data(:,2)>10);
Sur_Atom=X;
Sur_Res=A.data(X,1);
Sur_Res=union(Sur_Res,[]);