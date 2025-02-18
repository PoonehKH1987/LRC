function[Prot] = Pro_Graph(Pro)
% Shortcut summary goes here
 XYZ=[Pro.Atom_X,Pro.Atom_Y,Pro.Atom_Z]; 
%dt=delaunay(XYZ); %Delaunay Triangulation
dt=delaunay3(Pro.Atom_X,Pro.Atom_Y,Pro.Atom_Z);
e=squareform(pdist(XYZ,'euclidean')); %calculating the Distince between each 2 atoms
r=1;
while(r~=size(dt,1))
 Min=6;
 mt1=zeros(4,4);
 m1=dt(r,:);
  for i=1:4
    j=i;
    for k=j+1:4
      mt1(i,k)=e(m1(i),m1(k));
      if(Min>mt1(i,k))
         Min=mt1(i,k);
         mindis(r,1)=m1(i);
         resSeq(r,1)=Pro.Atom_resSeq1(m1(i),1);
         mindis(r,2)=m1(k);   
         Dis_res1(r,1)=mt1(i,k);
         resSeq(r,2)=Pro.Atom_resSeq2(1,m1(k));
         resName(r,1)=Pro.Atom_resName1(m1(i),1);
        resName(r,2)=Pro.Atom_resName2(1,m1(k));
      end     
   end
  end
 r=r+1;
end

k=1;
for i=1:size(mindis,1)
    if mindis(i,1)~=0
        Mindis(k,:)=mindis(i,:);
        ResSeq(k,:)=resSeq(i,:);
        ResName(k,:)=resName(i,:);
        k=k+1;
    end
end
for j=1:size(Mindis,1)
    if ResSeq(j,1)==ResSeq(j,2) 
        Mindis(j,:)=0;
    end
end

for j=1:size(Mindis,1)
    res1=ResSeq(j,1);
    res2=ResSeq(j,2);
    for k=j+1:size(Mindis,1)
        if (res1==ResSeq(k,1) & res2==ResSeq(k,2)) | (res1==ResSeq(k,2) & res2==ResSeq(k,1))
          Mindis(k,:)=0;
    end
    end
end

k=1;
for i=1:size(Mindis,1)
    if Mindis(i,1)~=0
        Mindis_M(k,:)=Mindis(i,:);
        ResSeq_M(k,:)=ResSeq(i,:);
        ResName_M(k,:)=ResName(i,:);
        k=k+1;
    end
end
M1=zeros(size(Pro.SurfaceRes_HorizonTal,1),size(Pro.SurfaceRes_HorizonTal,1));
for i=1:size(M1,1)
    for j=1:size(M1,1)
        for k=1:size(ResSeq_M,1)
           if (Pro.SurfaceRes_HorizonTal(i,1)==ResSeq_M(k,1) & Pro.SurfaceRes_VerTal(1,j)==ResSeq_M(k,2)) | (Pro.SurfaceRes_HorizonTal(i,1)==ResSeq_M(k,2) & Pro.SurfaceRes_VerTal(1,j)==ResSeq_M(k,1))
             M1(i,j)=1;
             M1(j,i)=1;
           end
        end
        
    end
end
Prot=Pro;
Prot.Matrix=M1;
%xlswrite(strcat('E:\My Thesis\',ProName{1,1,'.xlsx'),Prot.Matrix,1);
%return Pro.Matrix 
end


