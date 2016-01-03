function [X,Y,Points,NumberOfPointsInCell]=ReadTess(filename)
% Read the HDF5 data
X=h5read(filename,'/geometry/x_coordinate');
Y=h5read(filename,'/geometry/y_coordinate');
NumberOfCells=length(X);
Vertx=h5read(filename,'/geometry/x_vertices');
Verty=h5read(filename,'/geometry/y_vertices');
nVert=h5read(filename,'/geometry/n_vertices');
maxfaces=max(nVert);
TotalVertices=1;
NumberOfPointsInCell=zeros(NumberOfCells,1);
Points=zeros(NumberOfCells,maxfaces,2);
for i=1:NumberOfCells
    n=nVert(i);
    Points(i,1:n,1)=Vertx(TotalVertices:TotalVertices+n-1);
    Points(i,1:n,2)=Verty(TotalVertices:TotalVertices+n-1);
    NumberOfPointsInCell(i)=n;
    TotalVertices=TotalVertices+n;
end
