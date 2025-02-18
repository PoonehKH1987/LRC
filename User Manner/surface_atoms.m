function[Sur_Res,Sur_Atom]=surface_atoms(ASA)
A=Importdata(char(ASA));
X=find(A.data(:,2)>10);
Sur_Atom=A.data(X,1);
Sur_Res=union(Sur_Res,[]);