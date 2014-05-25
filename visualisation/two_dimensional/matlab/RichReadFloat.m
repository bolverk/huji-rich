function [X,Y,Pressure,Density,xVelocity,yVelocity,Points,time,Tracers,NumberOfPointsInCell]=RichReadFloat(filename,ShouldPlot,WhatToPlot,LogScale)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Matlab script to read RICH binary files in float format
%
%   Inputs  : filename - The location of the binary file 
%             ShouldPlot - 2: indicates to the current figure, 1: indicates
%             that a new figure is wanted 0: no plot
%             WhatToPlot - 1: density, 2: pressure, 3: pressure/density, 4:
%             Entropy (assuming gamma=5/3), 5: Tracer (assuming
%             tracerindex=1)
%             LogScale - 1 Plot Log scaled, 0 linear scale
%   Output  : X - The x cordinate mesh points
%             Y - The y cordinate mesh points
%             Pressure - The pressure
%             Density - The density
%             xVelocity - The x component of the velocity
%             yVelocity - The y component of the velocity
%             Points - The vertices of the Voronoi cell
%             time - The time of the simulation
%             Tracers - The scalar tracers
%             NumberOfPointsInCell - The number of vertices in each Voronoi
%             cell
%   Example : [X,Y,Pressure,Density,Points,xVelocity,yVelocity,time,Tracers,NumberOfPointsInCell]=RichReadFloat("output.bin",1,1,1)
%               Plots a log scaled density plot from the file "output.bin"
%             [X,Y,Pressure,Density,Points,xVelocity,yVelocity,time,Tracers,NumberOfPointsInCell]=RichReadFloat("output.bin",0)
%               Only reads the data from "output.bin"
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (nargin==2),
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

draw=WhatToPlot;
Log=LogScale;
fid=fopen(filename);
N=fread(fid,1,'int');
data=fread(fid,N*2,'float');
X=reshape(data,2,N);
Y=X(2,:)';
X=X(1,:)';
Pressure=fread(fid,N,'float');
Density=fread(fid,N,'float');
xVelocity=fread(fid,N,'float');
yVelocity=fread(fid,N,'float');
Temperature=Pressure./Density;
maxfaces=14;

if(ShouldPlot==1||ShouldPlot==2)
    Vertices=zeros(maxfaces*N,2);
    Faces=zeros(maxfaces,N);
end
NumberOfPointsInCell=zeros(N,1);
Points=zeros(N,maxfaces,2);
for i=1:N
    n=fread(fid,1,'int');
    data=fread(fid,2*n,'float');
    data=reshape(data,2,n);
    data=data';
    if(ShouldPlot==1||ShouldPlot==2)
        Vertices((i-1)*maxfaces+1:(i-1)*maxfaces+n,:)=data;
        Faces(1:n,i)=((i-1)*maxfaces+1):((i-1)*maxfaces+n);
        Faces(n+1:maxfaces,i)=NaN;
    end
    Points(i,1:n,:)=data;
    NumberOfPointsInCell(i)=n;
end
time=fread(fid,1,'float');
TracerLength=fread(fid,1,'int');
TracerDim=fread(fid,1,'int');
Tracers=fread(fid,TracerLength*TracerDim,'float');
Tracers=reshape(Tracers,TracerDim,TracerLength);
fclose(fid);
if(ShouldPlot==1||ShouldPlot==2)
    if(ShouldPlot==2)
        f1=gcf;
    else
        f1=figure;
    end
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
    axis equal;
end
