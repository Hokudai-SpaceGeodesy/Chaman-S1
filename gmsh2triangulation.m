function TR = gmsh2triangulation( filename )
%GMSH2TRIANGULATION Read Gmsh-based mesh file (xxx.msh) and generate Triangulation class element
%   to be used in greenlinear_tde2 and smoothing_tde.
%   Mesh file was generated Gmsh ver 4.
%   Originally, gmsh2trirep was based on Gmsh ver 2 and obsolate TriRep class 
%   TR.Points           <--- TR.X (TriRep style)
%   TR.ConnectivityList <--- TR.Triangulation (TriRep style)
%   copyright (c) Masato Furuya  May 2020

fid = fopen(filename);

meshformat0 = textscan(fid, '%s',1);
meshformat1 = textscan(fid, '%f %d %d',1);
meshformat2 = textscan(fid, '%s',1);

Entities = textscan(fid, '%s',1);
Nen = textscan(fid, '%d %d %d %d',1);
Np = Nen{:,1};Nl = Nen{:,2};Nsuf = Nen{:,3};Nloop = Nen{:,4};
for i=1:Np
    xp = textscan(fid,' %d %f %f %f %d',1);
end

for i=1:Nl
    yp = textscan(fid,' %d %f %f %f %f %f %f %d %d %d %d',1);
end

for i=1:Nsuf
 Sp = textscan(fid,' %d %f %f %f %f %f %f %d %d %d %d %d %d %d',1);
end

EndEntities=textscan(fid, '%s',1);
%fgets(fid)
Nodes=textscan(fid, '%s',1);
%fgets(fid)

Nn = textscan(fid, '%d %d %d %d',1);
Nblk=Nn{:,1};
Nnodes=Nn{:,2}

X=[];P=[];
for i=1:Nblk
 Blk=textscan(fid,'%d %d %d %d',1);
 nn=Blk{:,4};
 for i=1:nn
     newnode=textscan(fid,'%d',1);
     X=[X;newnode{:,1}];
 end
 for i=1:nn
     nodepos=textscan(fid,'%f %f %f',1);
     P=[P;nodepos{:,1} nodepos{:,2} nodepos{:,3}];
 end
end

EndNodes=textscan(fid, '%s',1); %fgets(fid);
Elements=textscan(fid, '%s',1); %fgets(fid);

T=[];
Nelem = textscan(fid, '%d %d %d %d',1);
nelem = Nelem{:,2}
Selem = textscan(fid, '%d %d %d %d',1);
for i=1:nelem
    Connect=textscan(fid,'%f %f %f %f',1);
    T=[T;Connect{:,2} Connect{:,3} Connect{:,4}];
end

TR=triangulation(T,P);

fclose(fid);

end

