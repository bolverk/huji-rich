function [X,Y,Pressure,Density,Vx,Vy,Points,time,Tracers,NumberOfPointsInCell]=read_hdf(filename,ShouldPlot,WhatToPlot,LogScale)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Matlab script to read RICH binary files in float/double format
%
%   Inputs  : filename - The location of the binary file 
%             ShouldPlot - 1 indicates that a plot is wanted 0 otherwise
%             WhatToPlot - 1: density, 2: pressure, 3: pressure/density, 4:
%             Entropy (assuming gamma=5/3), 5: Tracer (assuming
%             tracerindex=1)
%             LogScale - 1 Plot Log scaled, 0 linear scale
%   Output  : X - The mesh points
%             Pressure - The pressure
%             Density - The density
%             Vx - The x component of the velocity
%             Vy - The y component of the velocity
%             Points - The vertices of the Voronoi cell
%             time - The time of the simulation
%             Tracers - The scalar tracers
%             NumberOfPointsInCell - The number of vertices in each Voronoi
%             cell
%   Example : [X,Pressure,Density,Points,xVelocity,yVelocity,time,Tracers,NumberOfPointsInCell]=read_hdf("output.bin",1,1,1)
%               Plots a log scaled density plot from the file "output.bin"
%             [X,Pressure,Density,Points,xVelocity,yVelocity,time,Tracers,NumberOfPointsInCell]=read_hdf("output.bin",0)
%               Only reads the data from "output.bin"
%             [X,Pressure,Density,Points,xVelocity,yVelocity,time,Tracers,NumberOfPointsInCell]=read_hdf("output.bin")
%               Only reads the data from "output.bin"
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(nargin==1),
ShouldPlot=0;
WhatToPlot=1;
LogScale=0;
elseif (nargin==2),
  ShouldPlot=0;
  WhatToPlot=1;
  LogScale=0;
elseif (nargin==3),
   LogScale=0;
elseif (nargin==4),
   % do nothing
else
   error('Illigal number of input arguments');
end

% Read the HDF5 data
Density=h5read(filename,'/density');
Pressure=h5read(filename,'/pressure');
X=h5read(filename,'/x_coordinate');
Y=h5read(filename,'/y_coordinate');
Vx=h5read(filename,'/x_velocity');
Vy=h5read(filename,'/y_velocity');
time=h5read(filename,'/time');
NumberOfTracers=h5read(filename,'/Number of tracers');
NumberOfCells=length(Density);
Vertx=h5read(filename,'/x position of vertices');
Verty=h5read(filename,'/y position of vertices');
nVert=h5read(filename,'/Number of vertices in cell');
Tracers=zeros(NumberOfCells,NumberOfTracers);

for i=1:NumberOfTracers
    Tracers(:,i)=h5read(filename,sprintf('/Tracer number %d',i));
end

draw=WhatToPlot;
Log=LogScale;
Temperature=Pressure./Density;
maxfaces=max(nVert);
if(maxfaces>14)
    display('Warning, max number of faces exceeds 14!!')
end

if(ShouldPlot==1)
    Vertices=zeros(maxfaces*NumberOfCells,2);
    Faces=zeros(maxfaces,NumberOfCells);
end
NumberOfPointsInCell=zeros(NumberOfCells,1);
Points=zeros(NumberOfCells,maxfaces,2);
TotalVertices=1;
for i=1:NumberOfCells
    n=nVert(i);
    if(ShouldPlot==1)
        Vertices((i-1)*maxfaces+1:(i-1)*maxfaces+n,1)=Vertx(TotalVertices:TotalVertices+n-1);
        Vertices((i-1)*maxfaces+1:(i-1)*maxfaces+n,2)=Verty(TotalVertices:TotalVertices+n-1);
        Faces(1:n,i)=((i-1)*maxfaces+1):((i-1)*maxfaces+n);
        Faces(n+1:maxfaces,i)=NaN;
    end
    Points(i,1:n,1)=Vertx(TotalVertices:TotalVertices+n-1);
    Points(i,1:n,2)=Verty(TotalVertices:TotalVertices+n-1);
    NumberOfPointsInCell(i)=n;
    TotalVertices=TotalVertices+n;
end

if(ShouldPlot==1)
    f1=figure;
    set(f1,'Units','normalized')
    set(f1, 'Position', [0.03 0.03 0.65 0.85])
    hold on;
    switch (draw)
        case 1
            if(Log==1)
                caxis([min(log10(Density)) max(log10(Density*1.01))]);
                patch('Faces',Faces','Vertices',Vertices,'FaceVertexCData',log10(Density),'FaceColor','flat','EdgeAlpha',0.05);
            else
                caxis([min(Density) max(Density)*1.01]);
                patch('Faces',Faces','Vertices',Vertices,'FaceVertexCData',Density,'FaceColor','flat','EdgeAlpha',0.05);
            end
        case 2
            if(Log==1)
                caxis([min(log10(Pressure)) max(log10(Pressure*1.01))]);
                patch('Faces',Faces','Vertices',Vertices,'FaceVertexCData',log10(Pressure),'FaceColor','flat','EdgeAlpha',0.05);
            else
                caxis([min(Pressure) max(Pressure)*1.01]);
                patch('Faces',Faces','Vertices',Vertices,'FaceVertexCData',Pressure,'FaceColor','flat','EdgeAlpha',0.05);
            end
        case 3
            if(Log==1)
                caxis([min(log10(Temperature)) max(log10(Temperature*1.01))]);
                patch('Faces',Faces','Vertices',Vertices,'FaceVertexCData',log10(Temperature),'FaceColor','flat','EdgeAlpha',0);
            else
                caxis([min(Temperature) max(Temperature)*1.01]);
                patch('Faces',Faces','Vertices',Vertices,'FaceVertexCData',Temperature,'FaceColor','flat','EdgeAlpha',0);
            end
        case 4
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Change here gamma as needed
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            gamma=5/3;
            Entropy=Pressure./Density.^gamma;
            if(Log==1)
                caxis([min(log10(Entropy)) max(log10(Entropy*1.01))]);
                patch('Faces',Faces','Vertices',Vertices,'FaceVertexCData',log10(Entropy),'FaceColor','flat','EdgeAlpha',0);
            else
                caxis([min(Entropy) max(Entropy)*1.01]);
                patch('Faces',Faces','Vertices',Vertices,'FaceVertexCData',Entropy,'FaceColor','flat','EdgeAlpha',0);
            end
        case 5
            tracerindex=1;
            if(Log==1)
                caxis([min(log10(Tracers(tracerindex,:))) max(log10(Tracers(tracerindex,:)*1.01))]);
                patch('Faces',Faces','Vertices',Vertices,'FaceVertexCData',log10(Tracers(tracerindex,:)'),'FaceColor','flat','EdgeAlpha',0);
            else
                caxis([min((Tracers(tracerindex,:))) max((Tracers(tracerindex,:)*1.01))]);
                patch('Faces',Faces','Vertices',Vertices,'FaceVertexCData',Tracers(tracerindex,:)','FaceColor','flat','EdgeAlpha',0.05);
            end
    end
    colorbar;
end
