function [E]=LRC(ASA,Antigen,chainID,X)
%ASA=The name of output file of GetArea
%Antigen=The name of PDB file of antigen
%chainID=The chain of antigen
%X=the matrix n*7 contain the criterion file
[Sur_Res,Sur_Atom]=surface_residues(ASA);
Protein=Antigen_Surface_Graph4(Sur_Atom,Sur_Res,Antigen,chainID);
Un_W_graph=Protein.Matrix;
[W,W_graph]=weigthed_Gragh(Un_W_graph,X);
[Cluster]=MCL(W_graph,Protein.SurfaceRes_HorizonTal,30);
[E]=Filtering(Cluster,W);
fprintf('%i\n',E)
